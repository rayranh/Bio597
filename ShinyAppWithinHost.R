library(deSolve)
library(shiny) 
library(shinythemes)


ui <- fluidPage(
  theme = shinytheme("yeti"),
  navbarPage("Marek's disease Within-Host Model",
             sidebarLayout(
               sidebarPanel(
                 sliderInput(inputId = "beta", label = "beta", value = 0.02, min = 0.00, max = 0.5, step = 0.01),
                 
                 sliderInput(inputId = "beta_2", label = "beta_2", value = 0.003, min = 0.00, max = 0.5, step = 0.001),
                 
                 sliderInput(inputId = "nu_A", label = "activation by B cell", value = 0.03, min = 0, max = 0.1, step = 0.01),
                 
                 sliderInput(inputId = "nu_B", label = "activation by T cell", value = 0.02, min = 0, max = 0.1, step = 0.01),
                 
                 sliderInput(inputId = "B_cells", label = "number of B cells", value = 50, min = 0, max = 100, step = 1),
                 
                 sliderInput(inputId = "T_cells", label = "number of T cells", value = 50, min = 0, max = 100, step = 1)
               ),
               
               mainPanel(
                 plotOutput("plot")
               )
             )
  )
)

#trying to add something 


# Define the system of differential equations (ODEs)
sir_equations <- function(time, variables, parameters) {
  with(as.list(c(variables, parameters)), {
    dB <- -M * B_cells - beta * Cb * B_cells - beta_2 * Ct * B_cells + (g1 * (Cb + Ct) / (g2 + (Cb + Ct)))
    dCb <- M * B_cells + beta * Cb * B_cells + beta_2 * Ct * B_cells - alpha * Cb
    dT <- -nu_A * Cb * T_cells - nu_b * Ct * T_cells + (h1 * (Cb + Ct) / (h2 + (Cb + Ct)))
    dAt <- nu_A * Cb * T_cells + nu_b * Ct * T_cells - beta_2 * Ct * At - beta * Cb * At
    dLt <- theta * (beta_2 * Ct * At + beta * Cb * At) - mu * Lt
    dCt <- (1 - theta) * (beta_2 * Ct * At + beta * Cb * At) - alpha_2 * Ct
    dZ <- mu * Lt
    df <- -nu_F * Lt * f
    dIf <- nu_F * Lt * f
    return(list(c(dB, dCb, dT, dAt, dLt, dCt, dZ, df, dIf)))
  })
}

server <- function(input, output) {
  
  # Define the parameters as reactive values 
  #changed to reactive to add B and T cells 
  parameters_values <- reactive({
    c(
      M = 0.1
      , beta = input$beta          #contact rate with B cells (every 34 hours/ 1.4days) 
      , beta_2 = input$beta_2         #contact rate with T cells (every 333 hour/ 13 days) 
      #, beta_3 = beta_2              #contact rate with Ct cells for At 41 days activated T cells leaving 
      #, beta_4 = beta              #contact rate with Cb cells for At 41 days activated T cells leaving
      , nu_A =  input$nu_A               #Activation rate with T cells CD4+ (333 hours/ 13days)
      , nu_b = input$nu_B                #Activation rate with B cells (166 hours/ 41days)
      , nu_F = 1/200                #Infection rate of follicular cells (166 hours/ 41days)
      , mu =   1/168           #Rate of Tumor Cells (every 72 hours)
      , alpha = 1/3                #death rate of B cells (every 33 hours)
      , alpha_2 = 1/50        #death rate of T cells (every 48 hours)
      , theta = 0.7                 #population of activated T cells 
      , g1 = 1/15            #incoming B cells (every 15 hours)
      , g2 =0.001    
      , h1 = 0                      #incoming T cells / determined no incoming T cells 
      , h2 = 10
    )
  })
  
  # Define initial values for the model
  initial_values <- reactive({c( 
    B_cells = input$B_cells,  
    Cb = 0, 
    T_cells = input$T_cells,
    At = 0,
    Lt = 0,
    Ct = 0,
    Z = 0,
    f = 5, 
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
  output$plot <- renderPlot({
    sol <- sir_values_1()  # Get the updated solution
    
    # Convert the solution to a data frame
    sol_df <- as.data.frame(sol)
    
    # Plot the data
    with(sol_df, {
      plot(time, B_cells, col = "black", type = "l", ylim = c(0, 100), xlab = "Time (Hours)", 
           ylab = "Population density", main = "Marek's Model", cex.lab = 1.5, xlim = c(0, 200))
      lines(time, T_cells, col = "red")
      lines(time, Cb, col = "green")
      lines(time, At, col = "blue")
      lines(time, Lt, col = "purple")
      lines(time, Ct, col = "yellow")
      lines(time, Z, col = "lightblue")
    })
    
    # Adding the legend
    legend("topright", legend = c("B_cells", "T_cells", "Cb", "At", "Lt", "Ct", "Z"),
           col = c("black", "red", "green", "blue", "purple", "yellow", "lightblue"), 
           lty = 1, bty = "n")
  })
  
}

# Finally, run the Shiny app
shinyApp(ui = ui, server = server)
