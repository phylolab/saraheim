#------------------
#Circular bar plots
#------------------

#
library(ggplot2)
library(dplyr)
setwd("set/you/own/path/here")

generate_polar_graph=function(dataset){
  dataset=dataset[order(dataset$class),]
  dataset$FC=abs(dataset$FC)
  
  # Generate colors linked to the p_value
  colors=rainbow(length(unique(dataset$class)))
  color_list=c()
  color_loop=0
  iterator=1
  prev_class=""

  p_values=dataset$BH
  max_p=max(-log(p_values))
  min_p=min(-log(p_values))
  range=max_p-min_p

  classes=dataset$class

  for(p_val in p_values){
    curr_class=classes[iterator]

    if(curr_class!=prev_class){ color_loop=color_loop+1; new_color=col2rgb(colors[color_loop]) }
    fract=((max_p+log(p_val))/range)*(1-0.1)+0.1 #  change the log p to a fraction between 0.1 and 1
    color_list=c(color_list,rgb(new_color[1]/255,new_color[2]/255,new_color[3]/255,fract))
    prev_class=curr_class
    iterator=iterator+1
  }
  dataset$color<-color_list

  #   Set a number of 'empty bar' to add at the end of each group
  # MODIFICATION
  # We do not put empty bars for filtered classes 
  # instead of:
  # to_add <- data.frame( matrix(NA, empty_bar*nlevels(dataset$class), ncol(dataset)))
  #to_add$class <- rep(levels(as.factor(dataset$class)), each=empty_bar)
  # we do:
  # to_add <- data.frame( matrix(NA, empty_bar*length(unique(dataset$class)), ncol(dataset)))
  # to_add$class <- rep(unique(as.factor(dataset$class)), each=empty_bar)
  # This solves the problem of the numbre of colors below
  
  empty_bar <- 3
  to_add <- data.frame( matrix(NA, empty_bar*length(unique(dataset$class)), ncol(dataset)) )
  colnames(to_add) <- colnames(dataset)
  to_add$class <- rep(unique(as.factor(dataset$class)), each=empty_bar)
  dataset <- rbind(dataset, to_add)
  dataset <- dataset %>% arrange(class)
  dataset$id <- seq(1, nrow(dataset))

  #   Get the name and the y position of each label
  label_data <- dataset
  number_of_bar <- nrow(label_data)
  angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar
  label_data$hjust <- ifelse( angle < -90, 1, 0)
  label_data$angle <- ifelse(angle < -90, angle+180, angle)

  #   Making segments representing the groups
  base_data <- dataset %>% 
    group_by(class) %>% 
    summarize(start=(min(id)-0.2), end=(max(id)+0.2) - empty_bar) %>% 
    rowwise() %>% 
    mutate(title=mean(c(start, end)))
  base_data$label=unique(label_data$class)

  #   Make the y_axis
  FC_vals=c(na.omit(dataset$FC))
  y_ranges=c(0:ceiling(max(FC_vals)))
  y_axis_df=data.frame(y_vals=y_ranges)
  y_axis_df$x_vals=rep(0,nrow(y_axis_df))
  y_angle=rep(360/nrow(dataset),nrow(y_axis_df))
  y_axis_df$angle=y_angle
  title_y=data.frame(x_vals=c(max(dataset$id)-1),y_vals=c(max(y_axis_df$y_vals)/2),labels=c("Fold Change"),angle=angle[max(dataset$id)])

  #   Make a ggplot graph
  p1 <- ggplot()+
  geom_bar(data=dataset, aes(group=class, fill=as.factor(id), y=FC, x=as.factor(id)), stat="identity") +
  scale_fill_manual(values=color_list) +
  coord_polar()+
  ylim(-5, (0.6+ceiling(max(na.omit(dataset$FC))))) +
    theme(
      #panel.background= element_rect(fill = "#DFFDFF"),
      panel.background= element_rect(fill = "#FFFFFF"),
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid = element_blank(),
      legend.position="none",
    )+

  #   Name of each bar
  geom_text(data=label_data, aes(x=id, y=FC+0.5, label=LipidAnn, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=2.5, angle= label_data$angle, inherit.aes = FALSE )+

  #   Colored segments
  geom_segment(data=base_data, aes(x = start, y = -0.5, xend = end, yend = -0.5), colour = colors, alpha=0.8, size=0.6 , inherit.aes = FALSE )+

  #   Y_axis
  geom_text(data=title_y, aes(x = x_vals, y = y_vals, label=labels), colour = "black", size=4, angle=title_y$angle, fontface="bold", inherit.aes = FALSE)+
  geom_segment(data=data.frame(),aes(x = 0, y = 0, xend = 0, yend =ceiling(max(na.omit(dataset$FC)))), color="black", alpha=0.8, size=0.6 , inherit.aes = FALSE )+
  geom_segment(data=data.frame(vals=c(0:ceiling(max(na.omit(dataset$FC))))),aes(x = -0.3, y = vals, xend = 0.3, yend =vals), color="black", alpha=0.8, size=0.6 , inherit.aes = FALSE )+
  geom_text(data=y_axis_df, aes(x = max(dataset$id), y = y_vals, label=y_vals), colour = "black", alpha=0.8, size=4, inherit.aes = FALSE)

  p2<-ggplot()+
    theme(
      panel.background= element_rect(fill = "#FFFFFF"),
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid = element_blank(),
      legend.position="none",
    )+
  geom_text(data=base_data, aes(y = c(1:nrow(base_data)), x = 0, label=label), colour = colors, alpha=0.8, size=4, fontface="bold", inherit.aes = FALSE)

  df_grey <- data.frame(sequence = seq(min_p,max_p,0.01))
  p3<-ggplot(df_grey, aes(x=sequence, y=sequence,color=sequence)) +
    geom_point(aes(colour = sequence)) +
      theme(
        panel.background= element_rect(fill = "#FFFFFF"),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        legend.position="right",
      )+
    scale_colour_gradient(low = "white", high = "black",name="-log(p-value)")

  return(list(p1,p2,p3))
}

#   Import the data

df_DFover=read.csv('W:/FBM/DBC/nsalamin/clownfish/D2c/sheim/data/Figures Linda/Tables for Linda/BH_07102021_p_adjusted_Lipido_merged_DFoverCF_ceramidechecksept.csv', sep=";")
df_GENover=read.csv('W:/FBM/DBC/nsalamin/clownfish/D2c/sheim/data/Figures Linda/Tables for Linda/BH_p_adjvalues_FC_Lipido_mzRT_GENoverSPEC2_07102021.csv', sep=";")

# tw ofirsts graphs
foldChange='Log2FC_DF.CF'
df_DFover$FC<-df_DFover[,foldChange]

filtered_df=df_DFover[df_DFover$FC > 2 & df_DFover$BH < 0.05 ,]
plots=generate_polar_graph(filtered_df)
pdf("DFover_sup_2.1_plot.pdf")
print(plots[1])
dev.off()
pdf("DFover_sup_2.1_names.pdf")
print(plots[2])
dev.off()
pdf("DFover_sup_2.1_p_val.pdf")
print(plots[3])
dev.off()

filtered_df=df_DFover[df_DFover$FC < (-2) & df_DFover$BH < 0.05 ,]
plots=generate_polar_graph(filtered_df)
pdf("DFover_inf_3.1_plot.pdf")
print(plots[1])
dev.off()
pdf("DFover_inf_3.1_names.pdf")
print(plots[2])
dev.off()
pdf("DFover_inf_3.1_p_val.pdf")
print(plots[3])
dev.off()

# Other graphs
foldChange='Log2FCGen.Spec'
df_GENover$FC<-df_GENover[,foldChange]

filtered_df=df_GENover[df_GENover$FC > 2 & df_GENover$BH < 0.05 ,]
plots=generate_polar_graph(filtered_df)
pdf("GENover_sup_2_plot.pdf")
print(plots[1])
dev.off()
pdf("GENover_sup_2_names.pdf")
print(plots[2])
dev.off()
pdf("GENover_sup_2_p_val.pdf")
print(plots[3])
dev.off()

filtered_df=df_GENover[df_GENover$FC < (-2) & df_GENover$BH < 0.05 ,]
plots=generate_polar_graph(filtered_df)
pdf("GENover_inf_2_plot.pdf")
print(plots[1])
dev.off()
pdf("GENover_inf_2_names.pdf")
print(plots[2])
dev.off()
pdf("GENover_inf_2_p_val.pdf")
print(plots[3])
dev.off()

