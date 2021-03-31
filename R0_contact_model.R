library(shiny)

set.seed(1103111)

df <- 1
t <- 1
cols <- c("black", "white",
          rep(c(rainbow(10), gray(level=seq(0.3, 0.8, 0.1)),"gold","brown"), 4))

# Function for generating coordinates
pts <- function(n){
  x <- sqrt(n)*1.1
  y <- n/x
  xy <- expand.grid(seq(0,1, length.out=ceiling(x)), seq(0,1,length.out=ceiling(y)))
  jit <- .2/sqrt(n)
  xy[,1] <- xy[,1]+ runif(nrow(xy), -jit, jit)
  xy[,2] <- xy[,2]+ runif(nrow(xy), -jit, jit)
  delete <- nrow(xy) - n 
  xy <- xy[-sample(nrow(xy), delete), ]
  xy <- data.frame(xy)
  names(xy) <- c("x","y")
  return(xy)
}

         
# shiny app
SI <- shinyApp(
  ui = (fluidPage(
    
    sidebarLayout(
      sidebarPanel(
                
       # store window height in dimensions  
       tags$head(tags$script(' 
                                var dimension = [600];
                                $(document).on("shiny:connected", function(e) {
                                    dimension[0] = window.innerHeight;
                                    Shiny.onInputChange("dimension", dimension);
                                });
                                $(window).resize(function(e) {
                                    dimension[0] = window.innerHeight;
                                    Shiny.onInputChange("dimension", dimension);
                                });
                            ')),
        
        p("Disease transmission "),
        sliderInput("N",
                    "N:",
                    min = 20,
                    max = 2000,
                    value = 200,
                    step=20),
        sliderInput("R0",
                    "R0:",
                    min = 0.75,
                    max = 10,
                    value = 1,
                    step=0.25),
        radioButtons("inf",
                     "Contact distribution",
                     c("constant"=0, "heterogeneous"=1),
                     inline=F, selected=0),
        actionButton("go", "Next time step", icon("arrow-alt-circle-right")),
        actionButton("clear", "Restart", icon("refresh")),
        hr(),
        sliderInput("vac",
                    "Vaccinated [%]:",
                    min = 0,
                    max = 100,
                    value = 0,
                    step=5),
        sliderInput("I0",
                    "Initial infected:",
                    min = 1,
                    max = 5,
                    value = 1,
                    step=1),
         width=2
        ),

      mainPanel(
        plotOutput("mainPlot"), width=10,
      )
    )
  )),
  
  server <- (function(input, output, session) {
    
    newPlot <- function(){
    n <- input$N
    pers <<- pts(n)
    pers$time <<- 0
    pers$time[as.logical(rbinom(n, 1, input$vac/100))] <<- -1
    pers$time[round(seq(1/(input$I0+1), 0.999, 1/(input$I0+1))*n)] <<- 1
    par(mar=c(0,0,0,0)+0.3)
    layout(matrix(c(1,1,1,1,2), nrow=1, byrow=T))
    plot(pers$x,pers$y, xlab="", ylab="", xaxt="n", yaxt="n", lwd=2, cex=8-sqrt(n/120),
         bg=cols[pers$time+2], pch=21)
    text(pers$x,pers$y, pers$time, cex=round(1+8/sqrt(input$N),1))
    t <<- 1

    }

    updatePlot <- function(){
       curInf <- which(pers$time == t)
        if(input$inf==0) newCont <-  sapply(curInf, function(x) 
                sample(setdiff(1:input$N, x), round(input$R0), replace = F))
       else {
          rcont <- rpois(length(curInf), input$R0)
          newCont <-  sapply(curInf, function(x) sample(setdiff(1:input$N, x), 
                             max(rcont), replace = F))
          if(!is.null(ncol(newCont)))
             for(i in 1:ncol(newCont))
                if(rcont[i]+1 <= max(rcont))
                   newCont[min((rcont[i]+1),max(rcont)):max(rcont) ,i] <- NA
       }
       # print(sum(!is.na(newCont))/length(curInf))
       t <<- t + 1
       newInf <- which(1:input$N %in% unique(as.numeric(newCont)) & pers$time == 0)
       pers$time[newInf] <<- t
       par(mar=c(0,0,0,0)+0.3)
       layout(matrix(c(1,1,1,1,2,1,1,1,1,3), nrow=2, byrow=T))
       plot(pers$x,pers$y, xlab="", ylab="", xaxt="n", yaxt="n", lwd=2, 
         cex=8-sqrt(input$N/120),
         bg=cols[pers$time+2], pch=21)
       text(pers$x,pers$y, pers$time, cex=round(1+8/sqrt(input$N),1))
       if(length(curInf) > 0){
          for(i in 1:length(curInf)){
           if(!is.null(dim(newCont)))
                 arrows(pers$x[curInf[i]], pers$y[curInf[i]],
                                   pers$x[newCont[,i]], pers$y[newCont[,i]],
                                   col="red", lwd=2, lty=1,
                                   angle=15, length = 0.3)
           else arrows(pers$x[curInf[i]], pers$y[curInf[i]],
                       pers$x[newCont[i]], pers$y[newCont[i]],
                       col="red", lwd=2, lty=1, angle=15, length = 0.3)
         }
       }
       It <- table(pers$time[pers$time > 0])
       ns <- sum(pers$time >= 0)
       par(mar=c(2,0,10,0)+0.3)
       plot(NA,NA, xlab="Time", ylab="", yaxt="n", xlim=c(1,max(17, t)), ylim=c(0,ns),
            cex.axis=2, cex.lab=2)
       abline(v=c(5,10,15), lty=3, col="gray")
       lines(1:max(pers$time), cumsum(It), col="red", lwd=3)
       lines(1:max(pers$time), ns-cumsum(It), col="darkgreen", lwd=3)
       text(max(pers$time), sum(It), sum(It), cex=3)
       text(max(pers$time), ns-sum(It), ns-sum(It), cex=3)
       par(mar=c(4,0,4,0)+0.3)
       barplot(table(pers$time[pers$time > 0]), col="red", xlim=c(1,max(17, t)*1.2),
               xaxt="n",ylim=c(0,ns/2.05),yaxt="n",xpd=NA,  xlab="Time", cex.lab=2)
       abline(h=round(seq(0,ns/2,length.out=3)), lty=3, col="gray")
       mtext(paste0(c(0,25,50),"%"), side=4, line=-2, cex=1.5,
           at=seq(0,ns/2,length.out=3))
    }

    observeEvent(input$go,{
          output$mainPlot <- renderPlot({
             updatePlot()
          }, height = ifelse(is.null(input$dimension[1]), 700, input$dimension[1]-6))
    })

    observeEvent(input$clear,{
          output$mainPlot <- renderPlot({
             newPlot()
          }, height = ifelse(is.null(input$dimension[1]), 700, input$dimension[1]-6))
    })
                             
    observeEvent(input$N,{
          output$mainPlot <- renderPlot({
             newPlot()
          }, height = ifelse(is.null(input$dimension[1]), 700, input$dimension[1]-6))
    })                         
    
    observeEvent(input$R0,{
      if(input$R0%%1!=0)
        updateRadioButtons(session, "inf", selected = 1)
      else
        updateRadioButtons(session, "inf", selected = 0)
    })

  })
  
)


runApp(SI, launch.browser=T)
