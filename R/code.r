library(gplots)
library(shiny)
library(shinyjs)
library(DT)
library(shinythemes)
library(reshape2)
library(ggplot2)
library(Gviz)
library(GenomicRanges)
library(rtracklayer)
library(tidyr)
library(visNetwork)

options(shiny.maxRequestSize=5000*1024^2)
options(encoding = "UTF-8")

ui <- fluidPage(
  useShinyjs(),  # Include shinyjs
  column(4,             br(),
         h1("Welcome to hicGraph",align="center"),
         h3("Please select your input data",align="center"),
         column(12,textInput("chr_num","Please input your chr_num:")),
         column(12,fileInput("file1", "Choose O/E matrix file (bin=80k)")),
         column(12,fileInput("file2", "Choose O/E matrix file (bin=20k)")),
         column(12,fileInput("file3", "Choose chipseq .bigwig file")),
         column(12,fileInput("file4", "Choose compartment bed file")),
         column(12,fileInput("file5", "Choose TAD bed file")),
         column(12,fileInput("file6", "Choose mcool files"))
         #column(12,textInput("file3","Input your chipseq .bigwig file dir："))
  ),
  column(8,br(),br(),column(6,uiOutput("compartmentUI",align="center")),
         column(6,uiOutput("TadUI",align="center")),
         column(width = 10,offset = 2,plotOutput("heatmap")),
         column(11,plotOutput("omicplot")),
         column(12,visNetworkOutput("netplot"))
  )
)

server <- function(input, output) {

  chr_num <- reactive({
    req(input$chr_num)
    input$chr_num
  })

  input_file1 <- reactive({
    inFile <- input$file1
    if (is.null(inFile)) {
      return("")
    }
    temp <- read.csv(input$file1$datapath,header = F)
    temp
  })

  input_file2 <- reactive({
    inFile <- input$file2
    if (is.null(inFile)) {
      return("")
    }
    temp <- read.csv(input$file2$datapath,header=T)
    temp
  })

  input_file3 <- reactive({
    inFile <- input$file3
    if(is.null(inFile)) {
      return("")
    }
    temp <- import(input$file3$datapath)
    temp
  })

  input_file4 <- reactive({
    inFile <- input$file4
    if (is.null(inFile)) {
      return("")
    }
    temp <- read.table(input$file4$datapath,header=T)
    temp
  })

  input_file5 <- reactive({
    inFile <- input$file5
    if (is.null(inFile)) {
      return("")
    }
    temp <- read.csv(input$file5$datapath,header=T,sep = "\t")
    temp
  })

  input_file6 <- reactive({
    inFile <- input$file6
    if (is.null(inFile)) {
      return("")
    }
    temp <- read.csv(input$file6$datapath,header=T,sep = "\t")
    temp
  })




  data1 <- reactive({
    chr_num <- chr_num()
    #读取TAD位置信息
    all_tad <- input_file5()
    chr_mcool <- input_file6()
    #读入compartment文件
    all_compartment <- input_file4()
    #提取对应染色体的compartment信息
    if((all_compartment)!=""&&(chr_mcool!="")&&(all_tad!=""))
    {
      chr_compartment <- all_compartment[which(all_compartment$chrom==chr_num),]
      chr_compartment$com_name <- paste("com_",c(1:nrow(chr_compartment)),sep="")

      #提取不同染色体的tad信息
      chr_tad <- all_tad[which(all_tad$chrom == chr_num),]
      chr_tad[,2:4] <- as.data.frame(lapply(chr_tad[,2:4],as.numeric))
      chr_tad$tad_name <- paste("tad_",c(1:nrow(chr_tad)),sep="")
      #去除端粒的交互信息
      min <- which(chr_mcool$start1==chr_compartment$start[1])[1]
      max <- tail(which(chr_mcool$end1==chr_compartment$end[nrow(chr_compartment)]),n=1)
      chr_mcool <- chr_mcool[min:max,]
      chr_mcool <- chr_mcool[order(chr_mcool$start2),]
      row.names(chr_mcool) <- c(1:nrow(chr_mcool))
      min <- which(chr_mcool$start2==chr_compartment$start[1])[1]
      max <- tail(which(chr_mcool$end2==chr_compartment$end[nrow(chr_compartment)]),n=1)
      chr_mcool <- chr_mcool[min:max,]
      row.names(chr_mcool) <- c(1:nrow(chr_mcool))
      chr_mcool <- chr_mcool[order(chr_mcool$start1),]

      chrom1_com <- c()#用于存储mcool文件中每行chrom1所在compartment
      chrom2_com <- c()#用于存储mcool文件中每行chrom2所在compartment
      m <- 1
      #通过循环，将mcool中每条交联信息位置对应到compartment位置
      for(i in 1:nrow(chr_mcool)){
        tad1 <- which(chr_mcool$start1[i] >= chr_compartment$start & chr_mcool$end1[i] <= chr_compartment$end)
        tad2 <- which(chr_mcool$start2[i] >= chr_compartment$start & chr_mcool$end2[i] <= chr_compartment$end)
        chrom1_com[m] <- paste('com_',tad1,sep="")
        chrom2_com[m] <- paste('com_',tad2,sep="")
        m <- m+1
      }
      #生成links文件
      compartment_links <- data.frame(chrom1_com,chrom2_com,chr_mcool$count)
      compartment_links <- tidyr::unite(compartment_links,"com",chrom1_com,chrom2_com,sep="-")
      compartment_links <- aggregate(compartment_links[,"chr_mcool.count"],by=list(compartment_links[,"com"]),mean)
      compartment_links <- separate(compartment_links,Group.1,c("source","target"),sep="-",remove = T)
      com_count <- length(unique(chrom1_com))
      colnames(compartment_links) <- c("source","target","value")
      #过滤掉link数小于平均值的线
      compartment_links <- compartment_links[which(compartment_links$value>mean(compartment_links$value)),]
      #给compartment分组
      groups <- chr_compartment$com
      groups[which(chr_compartment$E1>0)] <- "compartment_A"
      groups[which(chr_compartment$E1<0)] <- "compartment_B"
      #计算nodes数
      #点的大小以边的count值的加和来计算
      value <- com_count
      for(i in 1:length(com_count)){
        value[i] <- sum(compartment_links[which(com_count[i]==compartment_links$source),]$value)
      }
      #生成nodes文件
      compartment_nodes <-data.frame(id=chr_compartment$com,
                                     label=chr_compartment$com,
                                     group=groups,
                                     title=paste(chr_compartment$start,chr_compartment$end,sep="_"),
                                     value=value)
      #在nodes中加入位置信息，用于与热图进行对照
      #当选择对应compartment后，绘制图形的同时提取nodes中的位置信息
      #在热图中找到大概的对应位置，重新绘制热图
      #下同（tad_nodes）
      compartment_nodes$start <- chr_compartment$start
      compartment_nodes$end <- chr_compartment$end
      #生成edges文件
      compartment_edges <- data.frame(from=compartment_links$source,to=compartment_links$target)

      list(compartment_nodes,compartment_edges,chr_tad,chr_mcool,chr_compartment)

    }
  })


  data2 <- reactive({
    chr_num <- chr_num()
    thedata <- data1()
    chr_tad <- thedata[[3]]
    chr_mcool <- thedata[[4]]
    chr_compartment <- thedata[[5]]
    m <- 1
    chr_tad$tad_com <- c(1:nrow(chr_tad))
    for(i in 1:nrow(chr_tad)){
      if(chr_tad$end[i] < chr_compartment$start[m]){
        next
      }
      else if(chr_tad$start[i] >= chr_compartment$start[m]){
        while(chr_tad$start[i] >= chr_compartment$end[m]){
          if(m<nrow(chr_compartment)){
            m <- m+1
          }
          else{
            break
          }
        }
        if(chr_tad$end[i] <= chr_compartment$end[m]){
          chr_tad$tad_com[i] <- chr_compartment$com[m]
        }
        else if(chr_tad$end[i] > chr_compartment$end[m] & m<nrow(chr_compartment)){
          chr_tad$tad_com[i] <- "2_com"
          if(m < nrow(chr_compartment))
          {
            m <- m+1
          }
          else{
            m <- m
          }
        }
      }
    }
    chr_tad_com <- chr_tad[which(chr_tad$tad_com==input$compartment_name),]
    if(is.na(chr_tad_com$chrom[1])){ #如果compartment中不含tad，打印以下语句
      print("No tad totally belong to this compartment!")
    }
    else{
      #提取tad内的compartment
      min_1 <- min(which(chr_mcool$start1 == chr_tad_com$start[1]))
      max_1 <- max(which(chr_mcool$end1 == chr_tad_com$end[nrow(chr_tad_com)]))
      chr_mcool_tad <- chr_mcool[min_1:max_1,]
      chr_mcool_tad <- chr_mcool_tad[order(chr_mcool_tad$start2),]
      min_2 <- min(which(chr_mcool_tad$start2 == chr_tad_com$start[1]))
      max_2 <- max(which(chr_mcool_tad$end2 == chr_tad_com$end[nrow(chr_tad_com)]))
      chr_mcool_tad <- chr_mcool_tad[min_2:max_2,]
      #统计tad间的交互关系
      chrom1_tad <- c()
      chrom2_tad <- c()
      m <- 1
      for(i in 1:nrow(chr_mcool_tad)){
        tad1 <- which(chr_mcool_tad$start1[i] >= chr_tad_com$start & chr_mcool_tad$end1[i] <= chr_tad_com$end)
        tad2 <- which(chr_mcool_tad$start2[i] >= chr_tad_com$start & chr_mcool_tad$end2[i] <= chr_tad_com$end)
        chrom1_tad[m] <- chr_tad_com$tad_name[tad1]
        chrom2_tad[m] <- chr_tad_com$tad_name[tad2]
        m <- m+1
      }
      #生成links
      tad_links <- data.frame(chrom1_tad,chrom2_tad,chr_mcool_tad$count)
      tad_links <- tidyr::unite(tad_links,"tad",chrom1_tad,chrom2_tad,sep="-")
      tad_links <- aggregate(tad_links[,"chr_mcool_tad.count"],by=list(tad_links[,"tad"]),mean)
      tad_links <- separate(tad_links,Group.1,c("source","target"),sep="-",remove = T)
      tad_count <- length(unique(chrom1_tad))
      colnames(tad_links) <- c("source","target","value")
      #点的大小以边的count值的加和来计算（效果不明显）
      value <- chr_tad_com$tad_name
      for(i in 1:length(chr_tad_com$tad_name)){
        value[i] <- sum(tad_links[which(chr_tad_com$tad_name[i]==tad_links$source),]$value)
      }
      #生成nodes文件
      tad_nodes <-data.frame(id=chr_tad_com$tad_name,
                             label=chr_tad_com$tad_name,
                             title=paste(chr_tad_com$start,chr_tad_com$end,sep="_"),
                             value=value)
      tad_nodes$start <- chr_tad_com$start
      tad_nodes$end <- chr_tad_com$end
      #生edge文件
      tad_edges <- data.frame(from=tad_links$source,to=tad_links$target)
      tad_title <- paste("network_node_tads_",chr_num,sep="")
      tad_title <- paste(tad_title,"_",sep="")
      tad_title <- paste(tad_title,input$compartment_name,sep="")

      list(tad_nodes,tad_edges,chr_mcool_tad,tad_title,chr_tad_com,chr_tad)
    }
  })


  data3 <- reactive({
    chr_num <- chr_num()
    thedata <-data2()
    chr_mcool_tad <- thedata[[3]]
    chr_tad_com  <- thedata[[5]]
    min_1 <- min(which(chr_mcool_tad$start1==chr_tad_com$start[which(chr_tad_com$tad_name==input$tad_name)]))
    max_1 <- max(which(chr_mcool_tad$end1==chr_tad_com$end[which(chr_tad_com$tad_name==input$tad_name)]))
    chr_mcool_tad_2 <-chr_mcool_tad[min_1:max_1,]
    chr_mcool_tad_2 <- chr_mcool_tad_2[order(chr_mcool_tad_2$start2),]
    min_2 <- min(which(chr_mcool_tad$start2==chr_tad_com$start[which(chr_tad_com$tad_name==input$tad_name)]))
    max_2 <- max(which(chr_mcool_tad$end2==chr_tad_com$end[which(chr_tad_com$tad_name==input$tad_name)]))
    chr_mcool_tad_2 <- chr_mcool_tad_2[min_2:max_2,]

    #统计fragement之间的交互信息
    fragment_links <- tidyr::unite(chr_mcool_tad_2,"fragment",c("start1","end1"),sep="_")
    fragment_links <- tidyr::unite(fragment_links,"fragment_2",c("start2","end2"),sep="_")
    fragment <- unique(fragment_links$fragment)
    fragment_name <- paste("fragment",c(1:length(fragment)))
    for(i in 1:length(fragment)){
      fragment_links$fragment_name[match(fragment[i],fragment_links$fragment)] <- fragment_name[i]
      fragment_links$fragment_2_name[match(fragment[i],fragment_links$fragment_2)] <- fragment_name[i]
    }
    #点的大小以边的count值的加和来计算
    fragment_value <- fragment
    for(i in 1:length(fragment)){
      fragment_value[i] <- sum(fragment_links[which(fragment_name[i]==fragment_links$fragment_name),]$count)
    }
    #生成fragment图的nodes和edges
    fragment_nodes <- data.frame(id=fragment_name,
                                 label=fragment_name,
                                 title=fragment,
                                 vlaue=fragment_value)
    fragment_edges <- data.frame(from=fragment_links$fragment_name,
                                 to=fragment_links$fragment_2_name,
                                 value=(fragment_links$count)/100)
    fragment_title <- paste("network_node_bins_",chr_num,sep="")
    fragment_title <- paste(fragment_title,"_",sep="")
    fragment_title <- paste(fragment_title,input$compartment_name,sep="")
    fragment_title <- paste(fragment_title,"_",sep="")
    fragment_title <- paste(fragment_title,input$tad_name,sep="")
    list(fragment_nodes,fragment_edges,fragment_title)
  })


  data4 <- reactive({
    if(input$compartment_name != "" && input$tad_name == ""){
      chr_num <- chr_num()
      oe_file_20k <- input_file2()
      thedata <- data1()
      compartment_nodes <- thedata[[1]]
      ############### genrate location information ##############
      bin <- 20000
      names_20k <- seq(from=0,by=bin,length.out = nrow(oe_file_20k))
      rownames(oe_file_20k) <- names_20k
      colnames(oe_file_20k) <- names_20k
      ############## Select compartment domains #############
      start_compartment <- compartment_nodes$start[which(compartment_nodes$id == input$compartment_name)]
      end_compartment <- compartment_nodes$end[which(compartment_nodes$id == input$compartment_name)]
      oe_file_compartment <- oe_file_20k[(which(as.numeric(rownames(oe_file_20k)) >= as.numeric( start_compartment))[1]-1):
                                           (which(as.numeric(rownames(oe_file_20k)) >= as.numeric(end_compartment))[1]+1),
                                         (which(as.numeric(rownames(oe_file_20k)) >= as.numeric(start_compartment))[1]-1):
                                           (which(as.numeric(rownames(oe_file_20k)) >= as.numeric(end_compartment))[1]+1)]
      start_heatmap <- rownames(oe_file_compartment)[1]
      end_heatmap <- rownames(oe_file_compartment)[nrow(oe_file_compartment)]
      list(start_heatmap,end_heatmap,oe_file_compartment)
    }
  })

  data5 <- reactive({
    if(input$compartment_name != "" && input$tad_name != ""){
      chr_num <- chr_num()
      oe_file_20k <- input_file2()
      thedata <- data2()
      tad_nodes <- thedata[[1]]
      ############### genrate location information ##############
      bin <- 20000
      names_20k <- seq(from=0,by=bin,length.out = nrow(oe_file_20k))
      rownames(oe_file_20k) <- names_20k
      colnames(oe_file_20k) <- names_20k
      ############## Select tad domains #############
      start_tad <- tad_nodes$start[which(tad_nodes$id == input$tad_name)]
      end_tad <- tad_nodes$end[which(tad_nodes$id == input$tad_name)]
      oe_file_tad <- oe_file_20k[(which(as.numeric(rownames(oe_file_20k)) >= as.numeric(start_tad))[1]-1):
                                   (which(as.numeric(rownames(oe_file_20k)) >= as.numeric(end_tad))[1]+1),
                                 (which(as.numeric(rownames(oe_file_20k)) >= as.numeric(start_tad))[1]-1):
                                   (which(as.numeric(rownames(oe_file_20k)) >= as.numeric(end_tad))[1]+1)]
      start_heatmap <- rownames(oe_file_tad)[1]
      end_heatmap <- rownames(oe_file_tad)[nrow(oe_file_tad)]
      list(start_heatmap,end_heatmap,oe_file_tad)
    }
  })

  data6 <- reactive({
    if(!is.null(input_file3())){
      chr_num <- chr_num()
      gr <- input_file3()
      genome <- 'hg19'
      ######### get chr name #######
      chr <- as.character(unique(seqnames(gr)))
      ######## get track message ########
      track <- DataTrack(range = gr,
                         genome = genome,
                         type = "histogram",
                         col.histogram = "#2167a4",
                         fill.histogram = "#2167a4",
                         name = " ",
                         showAxis = TRUE)
      list(track)
    }
  })


  output$compartmentUI <- renderUI({
    chr_num <- chr_num()
    all_tad <- input_file5()
    chr_mcool <- input_file6()
    #读入compartment文件
    all_compartment <- input_file4()
    #提取对应染色体的compartment信息
    if((all_compartment)!=""&&(chr_mcool!="")&&(all_tad!=""))
    {
      tempdata <- data1()
      compartment_nodes <- tempdata[[1]]
      selectInput("compartment_name","Select the compartment",choices = c("",unique(compartment_nodes$id)))
    }else{
      selectInput("compartment_name","Select the compartment",choices = c(""))
    }
  })

  output$TadUI <- renderUI({
    chr_num <- chr_num()
    #提取对应染色体的compartment信息
    if((input$compartment_name!=""))
    {
      tempdata <- data2()
      tad_nodes <- tempdata[[1]]
      selectInput("tad_name","Select the tad",choices = c("",unique(tad_nodes$id)))
    }else{
      selectInput("tad_name","Select the tad",choices = c(""))
    }
  })



  output$heatmap <- renderPlot({
    chr_num <- chr_num()
    if(!is.null(input_file1()) && !is.null(input_file2())) {
      oe_file_80k <- input_file1()
      if (input$compartment_name=="" && input$tad_name=="") {
        ############### genrate location information ##############
        bin <- 80000
        names_80k <- seq(from=0,by=bin,length.out = nrow(oe_file_80k))
        rownames(oe_file_80k) <- names_80k
        colnames(oe_file_80k) <- names_80k
        ############## Matrix #################
        oe_file_80k <- as.matrix(oe_file_80k)
        ################# ggplot2 #################
        options(scipen=200)
        melted_oe_file_80k <- melt(oe_file_80k)
        title_80k <- paste(chr_num,"_80k_heatmap",sep="")
        ggplot(data = melted_oe_file_80k,aes(x=Var1,y=Var2,fill=log(value+1)))+
          geom_raster()+
          ggtitle(title_80k)+
          scale_fill_gradient(low = 'white', high = 'red', #设置热图颜色
                              #breaks = scales::breaks_extended(n=101),
                              #na.value = 'red',
                              oob = scales::squish, #设定超出值的颜色
                              limits = c(0,1))+  #设定颜色范围
          scale_x_continuous(name = "",expand = c(0,0),position = "top")+ #去掉X轴标题，与X轴间隙
          scale_y_continuous(name = "",expand = c(0,0))+ #去掉Y轴标题，与Y轴间隙
          theme(panel.grid.major = element_blank(),#去掉背景色、网格线
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                axis.ticks.y = element_blank(), # 去除y轴刻度标签
                axis.text.y = element_blank(), #去除y轴刻度线
                legend.title = element_blank(),
                #axis.line = element_line(colour = "black"),#添加坐标轴
                #legend.position = "none",#去掉图例
                plot.title = element_text(hjust=0.5))#调整标题位置

      }else if(input$compartment_name!=""&&input$tad_name=="") {
        thedata <- data4()
        start_heatmap <- as.numeric(thedata[[1]])
        end_heatmap <- as.numeric(thedata[[2]])
        oe_file_compartment <- thedata[[3]]
        ############## Matrix #################
        oe_file_compartment <- as.matrix(oe_file_compartment)
        ############## ggplot2 ################
        options(scipen=200)
        melted_oe_file_compartment <- melt(oe_file_compartment)
        title_20k <- paste(chr_num,"_20k_heatmap",sep="")
        heat_loc <- paste(start_heatmap,end_heatmap,sep="-")
        title_20k <- paste(title_20k,heat_loc,sep="_")
        ggplot(data = melted_oe_file_compartment,aes(x=Var1,y=Var2,fill=log(value+1)))+
          geom_raster()+
          ggtitle(title_20k)+
          scale_fill_gradient(low = 'white', high = 'red', #设置热图颜色
                              oob = scales::squish, #设定超出值的颜色
                              limits = c(0,1))+  #设定颜色范围
          scale_x_continuous(name = "",expand = c(0,0),position = "top")+ #去掉X轴标题，与X轴间隙
          scale_y_continuous(name = "",expand = c(0,0))+ #去掉Y轴标题，与Y轴间隙
          theme(panel.grid.major = element_blank(),#去掉背景色、网格线
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                axis.ticks.y = element_blank(), # 去除y轴刻度标签
                axis.text.y = element_blank(), #去除y轴刻度线
                legend.title = element_blank(),
                plot.title = element_text(hjust=0.5))#调整标题位置
      }else if(input$compartment_name!="" && input$tad_name!=""){
        thedata <- data5()
        start_heatmap <- as.numeric(thedata[[1]])
        end_heatmap <- as.numeric(thedata[[2]])
        oe_file_tad <- thedata[[3]]
        ############## Matrix #################
        oe_file_tad <- as.matrix(oe_file_tad)
        ################# ggplot2 #################
        options(scipen=200)
        melted_oe_file_tad <- melt(oe_file_tad)
        title_20k <- paste(chr_num,"_20k_heatmap",sep="")
        heat_loc <- paste(start_heatmap,end_heatmap,sep="-")
        title_20k <- paste(title_20k,heat_loc,sep="_")
        ggplot(data = melted_oe_file_tad,aes(x=Var1,y=Var2,fill=log(value+1)))+
          geom_raster()+
          ggtitle(title_20k)+
          scale_fill_gradient(low = 'white', high = 'red', #设置热图颜色
                              oob = scales::squish, #设定超出值的颜色
                              limits = c(0,1))+  #设定颜色范围
          scale_x_continuous(name = "",expand = c(0,0),position = "top")+ #去掉X轴标题，与X轴间隙
          scale_y_continuous(name = "",expand = c(0,0))+ #去掉Y轴标题，与Y轴间隙
          theme(panel.grid.major = element_blank(),#去掉背景色、网格线
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                axis.ticks.y = element_blank(), # 去除y轴刻度标签
                axis.text.y = element_blank(), #去除y轴刻度线
                legend.title = element_blank(),
                plot.title = element_text(hjust=0.5))#调整标题位置
      }
    }
  })

  output$omicplot <- renderPlot({
    if(!is.null(input_file3()) && input$compartment_name == "" && input$tad_name == ""){
      thedata <- data6()
      track <- thedata[[1]]
      chr_num <- chr_num()
      ######## plot ########
      plotTracks(trackList = track,
                 chromosome = chr_num,
                 background.panel = "#f6f6f6",
                 background.title = "white",
                 col.title = "black",col.axis = "black",
                 rot.title = 0,cex.title = 0.5,margin=10,title.width = 1.75,
                 cex.axis = 1)
    }else if(input$compartment_name!=""&&input$tad_name==""){
      thedata <- data4()
      start_heatmap <- as.numeric(thedata[[1]])
      end_heatmap <- as.numeric(thedata[[2]])
      chr_num <- chr_num()
      thedata_2 <- data6()
      track <- thedata_2[[1]]
      ######## plot ########
      plotTracks(trackList = track,
                 chromosome = chr_num,
                 background.panel = "#f6f6f6",
                 from = start_heatmap, to =end_heatmap,
                 background.title = "white",
                 col.title = "black",col.axis = "black",
                 rot.title = 0,cex.title = 0.5,margin=10,title.width = 1.75,
                 cex.axis = 1)
    }else if(input$compartment_name!="" && input$tad_name!=""){
      thedata <- data5()
      start_heatmap <- as.numeric(thedata[[1]])
      end_heatmap <- as.numeric(thedata[[2]])
      chr_num <- chr_num()
      thedata_2 <- data6()
      track <- thedata_2[[1]]
      ######## plot ########
      plotTracks(trackList = track,
                 chromosome = chr_num,
                 background.panel = "#f6f6f6",
                 from = start_heatmap, to =end_heatmap,
                 background.title = "white",
                 col.title = "black",col.axis = "black",
                 rot.title = 0,cex.title = 0.5,margin=10,title.width = 1.75,
                 cex.axis = 1)
    }
  })

  output$netplot <- renderVisNetwork({
    chr_num <- chr_num()
    if (!is.null(data1())) {
      if (input$compartment_name==""&&input$tad_name=="") {
        thedata = data1()
        compartment_nodes <- thedata[[1]]
        compartment_edges <- thedata[[2]]
        compartment_title <- paste("network_node_compartments_",chr_num,sep = "")
        compartment_network <- visNetwork(compartment_nodes,compartment_edges,main=compartment_title)%>%
          visOptions(selectedBy = "group")%>%
          visLegend(useGroups = TRUE,width = 0.3,position = "right")%>%
          visOptions(highlightNearest = TRUE)%>%
          visInteraction(dragNodes = TRUE,## 可移动节点
                         dragView = TRUE, # 移动图
                         zoomView = TRUE,  # 缩放
                         navigationButtons = TRUE # 下方添加按钮
          )%>%
          visLayout(randomSeed = 4)
        compartment_network
      }else if (input$compartment_name!=""&&input$tad_name=="") {
        thedata = data2()
        tad_nodes <- thedata[[1]]
        tad_edges <- thedata[[2]]
        tad_title <- thedata[[4]]

        tad_network <- visNetwork(tad_nodes,tad_edges,main=tad_title)%>%
          visInteraction(dragNodes = TRUE,## 可移动节点
                         dragView = TRUE, # 移动图
                         zoomView = TRUE,  # 缩放
                         navigationButtons = TRUE # 下方添加按钮
          )%>%
          visIgraphLayout(layout = "layout_in_circle")
        tad_network
      }else {
        thedata = data3()
        fragment_nodes <- thedata[[1]]
        fragment_edges <- thedata[[2]]
        fragment_title <- thedata[[3]]

        network <- visNetwork(fragment_nodes,fragment_edges,main=fragment_title)%>%
          visInteraction(dragNodes = TRUE,## 可移动节点
                         dragView = TRUE, # 移动图
                         zoomView = TRUE,  # 缩放
                         navigationButtons = TRUE # 下方添加按钮
          )%>%
          visLayout(randomSeed = 4)
        network
      }
    }

  })

}

#' @import visNetwork
#' @import shiny
#' @export
run_app <- function(){
  shinyApp(ui, server)
}

