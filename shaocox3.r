#************************************Shao Cox******************************************;
#**************************************************************************************;
#** data.input: Input dataframe that must contain:                                    **;
#**      1) event.time:                             Event Time var (in days say)     **;
#**      2) censor.ind{1=Event, 0=Censored}:        Censoring Indicator              **;
#**      3) trt.ind{1=Immediate 0=Delayed}:         Treatment variable               **;
#**      4) switch.ind{1=Switched, 0=Not-switched}: Switching indicator called       **;
#**      5) switch.time:                            Switch time (in days say)        **;
#**************************************************************************************; 
Shao_Cox<-function(data.input)
{
	data.input=data.input[which(data.input$event.time!=0),]
	nrow=dim(data.input)[1]
	if(nrow<2) stop("Wrong dimention of input data")
	#data for imputation
	imput.data <- NULL

	for(i in 1:nrow)
	{	
			if(data.input$switch.ind[i]==1)
			{   
				#before switch
				temp=c(0,data.input$switch.time[i],0,data.input$trt.ind[i],0,0,0,0)
				#temp=c(0,data.input$switch.time[i],0,data.input$trt.ind[i],0,0)
				imput.data=rbind(imput.data,temp)
				#after switch time
				sw_trt=data.input$trt.ind[i]*data.input$switch.time[i]
				sw_ctrl=(1-data.input$trt.ind[i])*data.input$switch.time[i]
				temp=c(data.input$switch.time[i],data.input$event.time[i],data.input$censor.ind[i],1-data.input$trt.ind[i],sw_trt,sw_trt^2,sw_ctrl,sw_ctrl^2)
				# temp=c(data.input$switch.time[i],data.input$event.time[i],data.input$censor.ind[i],1-data.input$trt.ind[i],sw_ctrl,sw_ctrl^2)
				
				imput.data=rbind(imput.data,temp)
		
			}else{
				temp=c(0,data.input$event.time[i],data.input$censor.ind[i],data.input$trt.ind[i],0,0,0,0)
				#temp=c(0,data.input$event.time[i],data.input$censor.ind[i],data.input$trt.ind[i],0,0)
				imput.data=rbind(imput.data,temp)
			}
	}	
	
	colnames(imput.data) = c("start","stop","event", "trtcov", "sw0cov", "sw0cov2", "sw1cov", "sw1cov2")
	# colnames(imput.data) = c("start","stop","event", "trtcov", "sw_ctrl", "sw_ctrl2")
  imput.data=data.frame(imput.data,check.names=FALSE, row.names=NULL)
	return(coxph(Surv(start, stop, event) ~ trtcov + sw0cov + sw0cov2 + sw1cov + sw1cov2, data=imput.data))
	# return(coxph(Surv(start, stop, event) ~ trtcov + sw_ctrl + sw_ctrl2, data=imput.data))
}
###########################################################################################################################################
