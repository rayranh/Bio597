library(deSolve)
library(shiny) 
library(shinythemes)


### User Interface ### 

ui <- fluidPage(
  theme = shinytheme("yeti"),
  navbarPage("Marek's disease Within-Host Model",
             sidebarLayout(
               sidebarPanel(
                 sliderInput(inputId = "beta", label = "B cells cytolytically infected by CB", value = 0.01, min = 0.00, max = 0.1, step = 0.01),
                 
                 sliderInput(inputId = "beta_2", label = "T cells cytolytically infected by CT", value = 0.01, min = 0.00, max = 0.1, step = 0.01),  
                 
                 sliderInput(inputId = "mu_o", label = "From Lymphoid Organ to PBL", value = 0.05, min = 0, max = 1, step = 0.01),
                 
                 sliderInput(inputId = "mu_p", label = "From PBL to Lymphoid Organ ", value = 0.05, min = 0, max = 1, step = 0.01),
                 
                 sliderInput(inputId = "nu_A", label = "activation by B cell", value = 0.05, min = 0, max = 0.2, step = 0.01),
                 
                 sliderInput(inputId = "nu_B", label = "activation by T cell", value = 0.05, min = 0, max = 0.15, step = 0.01), 
                 
                 sliderInput(inputId = "alpha", label = "B cell Death", value = 0.3, min = 0, max = 1, step = 0.01),
                 
                 sliderInput(inputId = "alpha_2", label = "T cell Death", value = 0.01, min = 0, max = 0.1, step = 0.01),  
                 
                 sliderInput(inputId = "alpha_B", label = "B&T cell death in Blood", value = 0.1, min = 0, max = 0.3, step = 0.01), 
                 
                 sliderInput(inputId = "g1", label = "How fast B cells recruited", value = 50, min = 0, max = 200, step = 10), 
                 
                 sliderInput(inputId = "g2", label = "half-max, half of g1", value = 5000, min = 0, max = 20000, step = 10),
                 
                 sliderInput(inputId = "h1", label = "How fast T cells recruited", value = 50, min = 0, max = 200, step = 10),  
                 
                 sliderInput(inputId = "h2", label = "half-max, half of h1", value = 5000, min = 5000, max = 20000, step = 10), 
                 
                 sliderInput(inputId = "epsilon", label = "feather follicle infection", value = 0.03, min = 0, max = 0.5, step = 0.05), 
                 
                 sliderInput(inputId = "alpha_3", label = "death of infected FFE ", value = 0.02, min = 0, max = 0.5, step = 0.05), 
                 
                 
               ),
               
               mainPanel( 
                 width = 6,
                 fluidRow(
                   column(12, plotOutput("plot1"))   # Plot 1 in its own row
                 ),
                 fluidRow(
                   column(12, plotOutput("plot2"))   # Plot 2 in its own row
                 ),
                 fluidRow(
                   column(12, plotOutput("plot3"))   # Plot 3 in its own row 
                 ), 
                 fluidRow(
                   column(12, plotOutput("plot4"))   # Plot 4 in its own row
                 )
               )
             )
  )
)


# Define the system of differential equations (ODEs)
sir_equations <- function(time, variables, parameters) {
  with(as.list(c(variables, parameters)), { 
    
    ### Bursa ### 
    dB <- -M * B_cells - beta * Cb * B_cells - beta_2 * Ct * B_cells + (g1 * (Cb + Ct) / (g2 + (Cb + Ct))) - mu_o*B_cells + mu_p*B_cells
    dCb <- M * B_cells + beta * Cb * B_cells + beta_2 * Ct * B_cells - mu_o*Cb + mu_p*Cbb_cells - alpha*Cb
    dT <- -M * T_cells - nu_A * Cb * T_cells - nu_b * Ct * T_cells + (h1 * (Cb + Ct) / (h2 + (Cb + Ct))) - mu_o*T_cells + mu_p*Tb_cells 
    dAt <- M * T_cells + nu_A * Cb * T_cells + nu_b * Ct * T_cells - beta_2 * Ct * At - beta * Cb * At - 2*M*At - mu_o*At + mu_p*Atb_cells
    dLt <- theta * (beta_2 * Ct * At + beta * Cb * At) - mu_o*Lt + mu_p*Ltb_cells + M*At - gamma*Lt 
    dCt <- (1 - theta) * (beta_2 * Ct * At + beta * Cb * At) - alpha_2 * Ct + M*At - mu_o*Ct + mu_p*Ctb_cells + M*At 
    dZ <- gamma*Lt 
    
    ### entering blood ###   
    dBb <- mu_o*B_cells + mu_o*B_Th - 2 * mu_p*Bb_cells - alpha_B*Bb_cells   
    dCbb <- mu_o*Cb + mu_o*Cb_Th - 2 * mu_p*Cbb_cells - alpha_B*Cbb_cells   # removed mu_o B cells b/c lung dynamics not being modelled 
    dTb <- mu_o*T_cells + mu_o*T_Th - 2 * mu_p*Tb_cells - alpha_B*Tb_cells  
    dAtb <- mu_o*At + mu_o*At_Th  - 2 * mu_p*Atb_cells # no death for At cells 
    dLtb <- mu_o* Lt_Th + mu_o * Lt - 2 * mu_p*Ltb_cells  # population of Lt?  
    dCtb <- mu_o*Ct + mu_o*Ct_Th - 2*mu_p*Ctb_cells    
    
    ### Thymus ### 
    dB_Th <- -M * B_Th  - beta * Cb_Th* B_Th - beta_2 * Ct_Th * B_Th + (g1 * (Cb_Th + Ct_Th) / (g2 + (Cb_Th + Ct_Th))) - mu_o*B_Th + mu_p*Bb_cells 
    dCb_Th <- M * B_Th + beta * Cb_Th * B_Th + beta_2 * Ct_Th * B_Th - alpha * Cb_Th - mu_o*Cb_Th + mu_p*Cbb_cells 
    dT_Th <- - M * T_Th -nu_A * Cb_Th * T_Th - nu_b * Ct_Th * T_Th + (h1 * (Cb_Th + Ct_Th) / (h2 + (Cb_Th + Ct_Th))) - mu_o*T_Th + mu_p*Tb_cells 
    dAt_Th <- M * T_Th + nu_A * Cb_Th * T_Th + nu_b * Ct_Th * T_Th - beta_2 * Ct_Th * At_Th - beta*Cb_Th*At_Th - 2*M*At_Th - mu_o*At_Th + mu_p * Atb_cells 
    dLt_Th <- theta * (beta_2 * Ct_Th * At_Th + beta * Cb_Th * At_Th) - mu_o*Lt_Th + mu_p*Ltb_cells + M*At_Th - gamma*Lt_Th 
    dCt_Th <- (1 - theta) * (beta_2 * Ct_Th * At_Th + beta * Cb_Th * At_Th) - alpha_2 * Ct_Th + M*At_Th - mu_o*Ct_Th + mu_p*Ctb_cells     
    dZ_Th <- gamma * Lt_Th 
    
    ### FFE ###  
    df <- -Ltb_cells*f*epsilon
    dIf <- Ltb_cells*f*epsilon - alpha_3*If 
    
    return(list(c(dB, dCb, dT, dAt, dLt, dCt, dZ, dBb, dCbb, dTb, dAtb, dLtb, dCtb, dB_Th, dCb_Th, dT_Th, dAt_Th, dLt_Th, dCt_Th, dZ_Th, df, dIf)))
  })
}

#### Server #### 

server <- function(input, output) {
  
  # Define the parameters as reactive values 
  #changed to reactive to add B and T cells 
  parameters_values <- reactive({
    c(
      M = 0.5
      , beta = input$beta          #contact rate with B cells (every 34 hours/ 1.4days) 
      , beta_2 = input$beta_2         #contact rate with T cells (every 333 hour/ 13 days)  
      , mu_o = input$mu_o
      , mu_p = input$mu_p
      , nu_A =  input$nu_A               #Activation rate with T cells CD4+ (333 hours/ 13days)
      , nu_b = input$nu_B                #Activation rate with B cells (166 hours/ 41days)
      , alpha = input$alpha                #death rate of B cells (every 33 hours)
      , alpha_2 = input$alpha_2        #death rate of T cells (every 48 hours) 
      , alpha_B = input$alpha_B 
      , gamma = 1/30                 #rate of tumor formation 
      , theta = 0.8                 #population of activated T cells 
      , g1 = input$g1          #incoming B cells (every 15 hours)
      , g2 = input$g2    
      , h1 = input$h1                      #incoming T cells / determined no incoming T cells 
      , h2 = input$h2 
      , epsilon = input$epsilon                      #feather follicle infection rate  
      , alpha_3 = input$alpha_3
    )
  })
  
  # Define initial values for the model
  initial_values <- reactive({c(  
    ### bursa ### 
    B_cells = 1000,  
    Cb = 0, 
    T_cells =10,
    At = 0,
    Lt = 0,
    Ct = 0,  
    Z = 0, 
    ### blood ### 
    Bb_cells = 100, 
    Cbb_cells = 0, 
    Tb_cells = 100,
    Atb_cells = 0,  
    Ltb_cells = 0, 
    Ctb_cells = 0, 
    ### THYMUS ### 
    B_Th = 10,
    Cb_Th = 0,
    T_Th = 1000,
    At_Th = 0,
    Lt_Th = 0,
    Ct_Th = 0, 
    Z_Th = 0, 
    ### FFE ### 
    f = 100, 
    If = 0
  )}) 
  # Solve the ODE model using the parameters
  sir_values_1 <- reactive({
    ode(
      y = initial_values(),
      times = seq(0, 200, by = 0.1), #changing sequential time to 0.1
      func = sir_equations,
      parms = parameters_values()
    )
  })
  
  # Render the plot
  # Render the first plot
  output$plot1 <- renderPlot({
    sol <- sir_values_1()
    sol_df <- as.data.frame(sol)
    
    with(sol_df, {
      plot(time, B_cells, col = "black", type = "l", ylim = c(0, 1000), xlab = "Time (days)",
           ylab = "Population density", main = "Bursa", cex.lab = 1.5, xlim = c(0, 500/24))
      lines(time, T_cells, col = "red")
      lines(time, Cb, col = "green") 
      lines(time, At, col = "blue") 
      lines(time, Lt, col= "purple") 
      lines(time, Ct, col = "yellow") 
      lines(time, Z, col = "orange") 
      legend("topright", legend = c("B_cells", "T_cells", "Cb", "At", "Lt", "Ct", "Z"),
             col = c("black", "red", "green", "blue", "purple", "yellow", "orange"), lty = 1, bty = "n")
    })
  })
  
  # Render the second plot
  output$plot2 <- renderPlot({
    sol <- sir_values_1() # storing my output into sol 
    sol_df <- as.data.frame(sol) #converting sol to df 
    
    with(sol_df, {
      plot(time, Bb_cells, col = "black", type = "l", ylim = c(0, 1000), xlab = "Time (days)",
           ylab = "Population density", main = "Peripheral Blood", cex.lab = 1.5, xlim = c(0, 500/24))
      lines(time, Cbb_cells, col = "red")
      lines(time, Tb_cells, col = "green") 
      lines(time, Atb_cells, col = "blue") 
      lines(time, Ltb_cells, col = "purple") 
      lines(time, Ctb_cells, col = "yellow") 
      legend("topright", legend = c("Bb", "Cbb", "Tb", "Atb", "Ltb", "Ctb"),
             col = c("black", "red", "green", "blue", "purple", "yellow"), lty = 1, bty = "n")
    })
  }) 
  
  output$plot3 <- renderPlot({
    sol <- sir_values_1()
    sol_df <- as.data.frame(sol)
    
    with(sol_df, {
      plot(time, B_Th, col = "black", type = "l", ylim = c(0, 1000), xlab = "Time (days)",
           ylab = "Population density", main = "Thymus", cex.lab = 1.5, xlim = c(0, 500/24))
      lines(time, T_Th, col = "red")
      lines(time, Cb_Th, col = "green") 
      lines(time, At_Th, col = "blue") 
      lines(time, Lt_Th, col= "purple") 
      lines(time, Ct_Th, col = "yellow") 
      lines(time, Z_Th, col = "orange") 
      legend("topright", legend = c("B_cells", "T_cells", "Cb", "At", "Lt", "Ct", "Z"),
             col = c("black", "red", "green", "blue", "purple", "yellow", "orange"), lty = 1, bty = "n")
    })
  })
    output$plot4 <- renderPlot({
    sol <- sir_values_1()
    sol_df <- as.data.frame(sol)
    
    with(sol_df, {
      plot(time, f, col = "black", type = "l", ylim = c(0, 1000), xlab = "Time (days)",
           ylab = "Population density", main = "Feather Follicles", cex.lab = 1.5, xlim = c(0, 500/24))
      lines(time, If, col = "red")
      legend("topright", legend = c("Feather Follicles", "Infected Feather Follicles"),
             col = c("black", "red"), lty = 1, bty = "n")
    })
  })
  
}


shinyApp(ui = ui, server = server)
