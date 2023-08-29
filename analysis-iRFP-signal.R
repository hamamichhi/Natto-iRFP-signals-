#
# 
# This work is licensed under a Creative Commons Attribution 4.0 International License.
#

library(tidyverse)
library(rstan)
library(bayesplot)
library(rethinking)
library(magrittr)
library(lubridate)
library(stringr)

analysis.name<-"Natto" 
datetime<-format(Sys.time(),"%Y%m%d%H%M%S") 



# make folders
dir.create(path=file.path("figure"), showWarnings = FALSE, recursive = TRUE)
dir.create(path=file.path("output"), showWarnings = FALSE, recursive = TRUE)
dir.create(path=file.path("results","sample_file"), showWarnings = FALSE, recursive = TRUE)

# load dataset
groups<-c("HCD","HVK","LVK","Normal")

df<-read_csv(file.path("data","iRFP-signal.csv"))%>%
	gather(key, value, -Strain)%>%
	transmute(strain=factor(str_remove(Strain,'[0-9].'), groups), 
		age=as.numeric(key), weight=value)

# visualize raw data
theme_base<-theme_bw()+theme(panel.background=element_blank(),
                             strip.background=element_blank(),
                             legend.background = element_blank(),
                             legend.margin = margin(2, 2, 2, 2),
                             legend.key.height=unit(0.2, "cm"))
pdf(file=file.path("figure", sprintf("%s.raw.%s.pdf", analysis.name, datetime)), 
	onefile=TRUE, width=2.5, height=2, pointsize=8, 
	bg = "white")

p<-ggplot(df)+theme_base
p<-p+geom_point(aes(x=age, y=weight, color=strain),
	shape=16, size=0.5, alpha=0.8)
#p<-p+scale_color_manual(values=QRFP_colors)
print(p)

dev.off()

# prepare modeling
data<-list(N=nrow(df), MAX_W=10, X=floor(df$age), Y=df$weight, S=as.integer(df$strain))
stanmodel <- stan_model(file=file.path("stan", "model_bw.stan"))

####################################
# result and analysis
####################################

analysis.name<-"Natto-iRFP"
datetime<-format(Sys.time(),"%Y%m%d%H%M%S")
calc.stats<-function(fit, pars){
  stats.df<-tibble()
  for (par in pars){
    data<-rstan::extract(fit, pars=par)[[1]]
    temp.df<-enframe(HPDI(as.vector(data)))%>%
      bind_rows(enframe(quantile(as.vector(data))))%>%
      mutate(par=par)
    stats.df%<>%bind_rows(temp.df)}
  return(stats.df)
}

# result of modeling
mcmc.file<-sprintf("%s.csv", analysis.name)
fit <- sampling(stanmodel, data=data,
  iter=20000, warmup=5000, thin=1, seed=20230821, chain=1,
  sample_file=file.path("results","sample_file", mcmc.file))


# Save computed statistics
# output MCMC diagnostics
color_scheme_set("viridis")
post<- as.matrix(fit) # single chains 
np <- nuts_params(fit)

# trace
pdf(file=file.path("figure", sprintf("%s.trace.%s.pdf", analysis.name, datetime)), 
	onefile=TRUE, width=10, height=5, pointsize=8, 
	bg = "white")
for (i in seq(dimnames(post)$parameters)){
	p<-mcmc_trace(post, pars=dimnames(post)$parameters[i], np=np,
		np_style = trace_style_np(div_color = "black", div_size = 0.5),
		size=0.1)
	print(p)
}
dev.off()

#neff
color_scheme_set("brightblue") # see help("color_scheme_set")
ratios <- neff_ratio(fit)
pdf(file=file.path("figure", sprintf("%s.neff.%s.pdf", analysis.name, datetime)), 
	onefile=TRUE, width=5, height=10, pointsize=8, 
	bg = "white")
mcmc_neff(ratios, size = 2)
dev.off()

# calc post stats
stats.df<-calc.stats(fit, pars=c("sigma1", "sigma2", "lp__"))%>%
	spread(key=name, value=value)%>%
	write_csv(file.path("output", 
		sprintf("stats.%s.%s.csv", analysis.name, datetime)))

ts.stats.df<-tibble()

temp<-rstan::extract(fit, pars="B")$B
for (wk in 1:10){
	for (st in 1:4){
		data<-temp[,wk,st]		
		temp.df<-enframe(HPDI(as.vector(data)))%>%
	  	bind_rows(enframe(quantile(as.vector(data))))%>%
	  	mutate(par="Base", week=wk, strain=factor(st, levels=1:4, labels=groups))
	  ts.stats.df%<>%bind_rows(temp.df)
	}
}

temp<-rstan::extract(fit, pars="est_Y")$est_Y
for (wk in 1:10){
	for (st in 1:4){
		data<-temp[,wk,st]		
		temp.df<-enframe(HPDI(as.vector(data)))%>%
	  	bind_rows(enframe(quantile(as.vector(data))))%>%
	  	mutate(par="Est", week=wk, strain=factor(st,  levels=1:4, labels=groups))
	  ts.stats.df%<>%bind_rows(temp.df)
	}
}

ts.stats.df%<>%spread(key=name, value=value)

ts.stats.df%>%
	write_csv(file.path("output", 
		sprintf("ts.stats.%s.%s.csv", analysis.name, datetime)))

# calc posterior distribution of difference of BW at given ages
temp<-rstan::extract(fit, pars="B")$B
ts.diff.df<-tibble()
for (wk in 1:10){
	data<-temp[,wk,]
	ts.diff.df%<>%bind_rows(tibble(week=wk, type="HCD-HVK", diff=data[,1]-data[,2]))
	ts.diff.df%<>%bind_rows(tibble(week=wk, type="HCD-LVK", diff=data[,1]-data[,3]))
	ts.diff.df%<>%bind_rows(tibble(week=wk, type="HCD-Normal", diff=data[,1]-data[,4]))
	ts.diff.df%<>%bind_rows(tibble(week=wk, type="HVK-LVK", diff=data[,2]-data[,3]))
}
diff.types<-c("HCD-HVK", "HCD-LVK", "HCD-Normal","HVK-LVK")
ts.diff.df%<>%mutate(type=factor(type, diff.types))

ts.diff.stats.df<-tibble()
for (wk in 1:10){
	for (tp in diff.types){
		data<-subset(ts.diff.df, (week==wk)&(type==tp))$diff		
		temp.df<-enframe(HPDI(as.vector(data)))%>%
	  	bind_rows(enframe(quantile(as.vector(data))))%>%
	  	mutate(week=wk, type=factor(tp,  levels=diff.types))
	  ts.diff.stats.df%<>%bind_rows(temp.df)
	}
}
ts.diff.stats.df%<>%spread(name, value)

ts.diff.stats.df%>%
	write_csv(file.path("output", 
		sprintf("ts.diff.stats.%s.%s.csv", analysis.name, datetime)))

# plot timeseries
pdf(file=file.path("figure",
	sprintf("ts.ststats.%s.%s.pdf", analysis.name, datetime)), 
	onefile=TRUE, width=3, height=2, pointsize=8, 
	bg = "white")

p<-ggplot(subset(ts.stats.df, par=="Base"))+theme_base
p<-p+geom_line(aes(x=week, y=`50%`, group=strain, color=strain),
	size=0.2)
p<-p+geom_ribbon(aes(x=week, ymin=`|0.89`, ymax=`0.89|`, 
	group=strain, fill=strain), alpha=0.2)
p<-p+coord_cartesian(xlim=c(1,10))
#p<-p+scale_color_manual(values=QRFP_colors)
#p<-p+scale_fill_manual(values=QRFP_colors)
p<-p+labs(y="signal", x="Time (week)")
print(p)

dev.off()

# plot timeseries Graphs of two groups
#for (lb in c("HVK", "LVK", "Normal")){
pdf(file=file.path("figure",
                   sprintf(paste(lb,"ts.ststats_two_data.%s.%s.pdf",sep="_"), analysis.name, datetime)), 
    onefile=TRUE, width=3, height=2, pointsize=8, 
    bg = "white")

p<-ggplot(subset(ts.stats.df, par=="Base" & strain %in% c("HCD",lb)))+theme_base
p<-p+geom_line(aes(x=week, y=`50%`, group=strain, color=strain),
               size=0.2)
p<-p+geom_ribbon(aes(x=week, ymin=`|0.89`, ymax=`0.89|`, 
                     group=strain, fill=strain), alpha=0.2)
p<-p+coord_cartesian(xlim=c(1,10))
p <- p + scale_x_continuous(breaks = seq(1, 10, by = 1))
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#p<-p+scale_color_manual(values=QRFP_colors)
#p<-p+scale_fill_manual(values=QRFP_colors)
p<-p+labs(y="signal", x="Time (week)")
print(p)

dev.off()
}



colors <- c("HCD" = "#2980B9", "HVK" = "#C0392B", "LVK" = "#27AE60", "Normal" = "#8E44AD")
for (lb in c("HVK", "LVK", "Normal")) {
  pdf(file=file.path("figure",
                     sprintf(paste(lb, "ts.ststats_two_data.%s.%s.pdf", sep = "_"), analysis.name, datetime)),
      onefile = TRUE, width = 3, height = 2, pointsize = 8,
      bg = "white")
  
  p <- ggplot(subset(ts.stats.df, par == "Base" & strain %in% c("HCD", lb))) + theme_base
  p <- p + geom_line(aes(x = week, y = `50%`, group = strain, color = strain), size = 0.2)
  p <- p + geom_ribbon(aes(x = week, ymin = `|0.89`, ymax = `0.89|`, group = strain, fill = strain), alpha = 0.4)
  p <- p + coord_cartesian(xlim = c(1, 10))
  p <- p + scale_x_continuous(breaks = seq(1, 10, by = 1)) # Weekly increments
  p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) # グリッドを消す
  p <- p + scale_color_manual(values = colors) # Apply darker color
  p <- p + scale_fill_manual(values = colors) # Apply darker color
  p <- p + labs(y = "signal", x = "Time (week)")
  print(p)
  
  dev.off()
}



# plot diff histogram
pdf(file=file.path("figure",
	sprintf("ts.diff.histogram.%s.%s.pdf", analysis.name, datetime)), 
	onefile=TRUE, width=1.6, height=2, pointsize=8, 
	bg = "white")

for (wk in 1:10){
	df<-subset(ts.diff.df, (week==wk))%>%
		select(diff, type)
	p<-ggplot(df)+theme_base
	p<-p+geom_histogram(aes(x=diff, y=..density..), bins=100, alpha=0.5)
	p<-p+labs(title=sprintf("week %s", wk))

	p<-p+geom_vline(xintercept=0, size=0.2, alpha=0.5, color="red")
	p<-p+geom_vline(data=subset(ts.diff.stats.df, (week==wk)),
		map=aes(xintercept=`|0.89`), size=0.2)
	p<-p+geom_vline(data=subset(ts.diff.stats.df, (week==wk)),
		map=aes(xintercept=`50%`), size=0.2, linetype="dashed")
	p<-p+geom_vline(data=subset(ts.diff.stats.df, (week==wk)),
		map=aes(xintercept=`0.89|`), size=0.2)

	p<-p+facet_grid(type~.)

	p<-p+labs(x="Difference", y="Density")

	print(p)
}

dev.off()

