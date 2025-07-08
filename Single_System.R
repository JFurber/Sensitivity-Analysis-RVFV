# Run simulation, taking x args where x = number of parameters analysed in SA

singleSim_OnlyCulex <- function(lambda_C, rho_max_C,
                           A_C, kappa_C) {
  
  #Parameters
  
  delta_t <<- 1 # In days, the time step
  N_ponds <<- 1 # How many ponds in each cell - not ever changed N_ponds from 1
  reduction_C <<- 0.5 # Based on Figure 8 in Reisen, W. K., Fang, Y., & Martinez, V. M. (2006). Effects of temperature on the transmission of west nile virus by Culex tarsalis (Diptera: Culicidae). Journal of medical entomology, 43(2), 309–17. Retrieved from http://www.ncbi.nlm.nih.gov/pubmed/16619616
  reduction_A <<- 0.58 # Based on Figure 8 in Reisen, W. K., Fang, Y., & Martinez, V. M. (2006). Effects of temperature on the transmission of west nile virus by Culex tarsalis (Diptera: Culicidae). Journal of medical entomology, 43(2), 309–17. Retrieved from http://www.ncbi.nlm.nih.gov/pubmed/16619616
  b_L1 <<- 1/(5*365) # birth rate for livestock
  mu_L1 <<- 1/(5*365) # death rate for livestock
  q_divided <<- 1E11 # parameter for the impact of the livestock on vector fecundity and gonotrophic cycles
  q_A <<- 0.007 # probability of transovarial transmission
  epsilon_L1 <<- 2/7# input by Sophie 1/2 - should it be 2/7?
  gamma_L1 <<- 1/30 # input by Sophie 1/3     - should it be 30?
  t_dep <<- 0.229 # average time for egg deposition = 0.229 days in lab ref(8) GL2018 supp
  
  days <- day_length
  
  # In order to make sure this file is consistent with master file - this does not change:
  
  if (modeloutput == "allMozadultsmeanlast28days") {
    
    # Load packages:
    library(deSolve)
    library(rootSolve)
    library(tidyverse)
    library(gridExtra)
    
    # Functions required:
    source("Functions/system_of_RVF_equations_DEN.R", chdir = TRUE)
    source("Functions/development_rate.R", chdir = TRUE)
    source("Functions/mortality_rate.R", chdir = TRUE)
    source("Functions/oviposition_rate.R", chdir = TRUE)
    source("Functions/Sharpe_DeMichele_parameters.R", chdir=TRUE)
    
    # SET PARAMETERS FROM ARGUMENT--------------------------------------------------
    
    assign("lambda_C", value = lambda_C, pos = 2, inherits = TRUE)
    assign("rho_max_C", value = rho_max_C, pos = 2, inherits = TRUE)
    assign("kappa_C", value = kappa_C, pos = 2, inherits = TRUE)
    assign("A_C", value = A_C, pos = 2, inherits = TRUE)
    
    # Temperature, Water Body Area (wbarea), Livestock Density - Set up
    
    database <- data.frame(daynum = 1:5000, tmean = constanttemp, wbarea1 = constantwb,
                           cellarea = cellarea, livestockda = livestocktotal)
    
    # propsearcharea represents the fraction of the total cell area that mosquitoes are likely to search in - they will search the whole cell
    database$propsearcharea <- 1
    # wbarea represents the portion of the water body area that mosquitoes are likely to interact with.
    database$wbarea <- database$wbarea1 * database$propsearcharea
    
    # RUN Culex SIMULATION ---------------------------------------------------------------
    
    # Only look at relevant data for the location(s) being studied - this is taken from the original code, but in essence it is the same spot each time
    database_thisloc <- database
    
    assign("database_thisloc", value = database_thisloc, pos = 2, inherits = TRUE)
    
    # Sequence from 0, 1, 2, 3,... until last day number:
    times <- seq(0, days, by = delta_t)
    
    # Additional functions:
    source("Functions/ponds.R", chdir = TRUE) 
    source("Functions/biting_rate.R", chdir = TRUE)
    
    # This function is not needed for this code as the temperature is set as a constant - this has been addressed in the system of RVF equations file.
    # source("Functions/temperature.R", chdir = TRUE) 
    
    # # State Variables:
        # State_RVF <- c(O_c = 20814563120,
        #                L_c = 217139002419,
        #                P_c = 63312191177,
        #                A1_c = 30712596078,
        #                F_c = 50471083343,
        #                A2_c = 632816460)

    # Used to obtain state variables:
    State_RVF <- c(O_c = 100,
                   L_c = 0,
                   P_c = 0,
                   A1_c = 0,
                   F_c = 0,
                   A2_c = 0)
    
    N_L1 <- livestocktotal
    
    # Parameters to feed into ode function:
    Pars_RVF <- c(rho_max_C, 
                  q_A, 
                  b_L1, 
                  mu_L1, 
                  epsilon_L1, 
                  gamma_L1, 
                  q_divided, 
                  as.numeric(N_L1), 
                  prop_find,
                  livestocktotal
                  )
    
    # Function for differential equations:
    out <- ode(y = State_RVF, 
               times = times, 
               func = system_of_RVF_equations_density_temp_Culex,
               parms = Pars_RVF,
               atol = 1e-1, rtol = 1e-1,
               maxsteps = 10000)
    
    out_RVF_limiting <<- as.data.frame(out)
    
    # Add results from this loc to vector:
    
    cell_ode_output <- out_RVF_limiting 
    
    # Take data from all locs and bind into single dataframe:
    data_output <- cell_ode_output
    
    #Since only looking at Culex, non-infectious, the column names have been changed to show this.
    colnames(data_output) <-
      c("Time",                                  
        "Eggs",                      
        "Larvae",                    
        "Pupae",                    
        "Nulliparous",       
        "Flyers",                    
        "Parous"   
)   
    
    # Shift columns to more intuitive order - this was :
    
    facetPlotData <- data_output %>%
      tidyr::pivot_longer(cols=c("Eggs",                       
                                 "Larvae",                     
                                 "Pupae",                      
                                 "Nulliparous",          
                                 "Flyers",                     
                                 "Parous"),                       
                          names_to = "Compartment",
                          values_to = "Value")
    
    facetPlotData$Compartment <- fct_inorder(facetPlotData$Compartment) # to keep them in order as above for plotting
    
    adults_cols <- c("Nulliparous",          
                     "Flyers",                     
                     "Parous")       
    
    toReturn <- data_output
    
    toReturn$allMozadults <-rowSums(toReturn[, adults_cols])
    
    # if chosen "one", then the return plot will be the state variable plot and show the population. 
    if (what == "one") {
      facetPlot <- ggplot(data = facetPlotData, aes(x = Time, y = Value)) +
        geom_line(size = 0.8) + 
        facet_wrap(~Compartment, scales = "free") + 
        theme_classic(base_size = 12) +
        labs(x = "Time (days)", y = "Population") + 
        theme(
          # plot.title = element_text(hjust = 0.5), strip.text = element_text(size = 7)
          plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),  # Centered and prominent title
          axis.title = element_text(size = 10),      # Axis titles are clear and legible
          axis.text = element_text(size = 8),       # Tick labels are a bit smaller
          strip.text = element_text(size = 7, face = "bold"),  # Facet labels are bold for clarity
          panel.grid.major = element_line(color = "gray90"),   # Optional: add light grid lines for reference
          panel.grid.minor = element_blank(),        # Remove minor grid lines to avoid clutter
          panel.spacing = unit(0.5, "lines")           # Adjust spacing between facets
        )
      
      facetPlot2 <- ggplot(data = facetPlotData, aes(x = Time, y = Value)) +
        geom_line(size = 0.8) +
        facet_grid(~Compartment, scales = "fixed") +
        theme_classic(base_size = 12) +
        labs(x = "Time (days)", y = "Population") +
        scale_y_continuous(
          breaks = seq(0, 2e11, by = 2.5e10),  # Manually specify evenly spaced breaks
          labels = scales::label_number(scale = 1, accuracy = 1e10)  # Format labels for clarity
        ) +
        theme(
          plot.title = element_text(hjust = 0.5, size = 8, face = "bold"),  # Centered and prominent title
          axis.title = element_text(size = 8),      # Axis titles are clear and legible
          axis.text = element_text(size = 6),       # Tick labels are a bit smaller
          strip.text = element_text(size = 5, face = "bold"),  # Facet labels are bold for clarity
          panel.grid.major = element_line(color = "gray90"),   # Optional: add light grid lines for reference
          panel.grid.minor = element_blank(),        # Remove minor grid lines to avoid clutter
          panel.spacing = unit(1, "lines")           # Adjust spacing between facets
        )
      facetPlot2
      
      ggsave(plot = facetPlot,
             filename = paste0("Output/", "Culex_SingleSim",  "_facetPlot.png"),
             width = 6, height = 4)
      ggsave(plot = facetPlot,
             filename = paste0("Output/", "Culex_SingleSim",  "_facetPlot.pdf"),
             width = 6, height = 3,device=pdf)
      
      ggsave(plot = facetPlot2,
             filename = paste0("Output/", "Culex_SingleSim",  "_facetPlot2.png"),
             width = 6, height = 3)
      ggsave(plot = facetPlot2,
             filename = paste0("Output/", "Culex_SingleSim",  "_facetPlot2.pdf"),
             width = 6, height = 2,device=pdf)
      
      write.csv(data_output, file = paste0("Output/", "Culex_SingleSim", ".csv"))
      
    }
    
    # If chosen "SA", i.e. sensitivity analysis, then the return is the sum of all adults in past year.
    if (what == "SA") {
      
      startday <- day_length-365
      mozMeanLast4weeks <- toReturn %>%
        dplyr::filter(Time > startday)
      return(sum(mozMeanLast4weeks$allMozadults))
    }
    
  }
  
  else {stop("Error - Model output in SA_SingleSimulation.R does not match mastercode")}
}

singleSim_OnlyAedes <- function(lambda_A, rho_max_A,
                                A_A, kappa_A) {
  
  #Parameters
  
  delta_t <<- 1 # In days, the time step
  N_ponds <<- 1 # How many ponds in each cell - not ever changed N_ponds from 1
  reduction_C <<- 0.5 # Based on Figure 8 in Reisen, W. K., Fang, Y., & Martinez, V. M. (2006). Effects of temperature on the transmission of west nile virus by Culex tarsalis (Diptera: Culicidae). Journal of medical entomology, 43(2), 309–17. Retrieved from http://www.ncbi.nlm.nih.gov/pubmed/16619616
  reduction_A <<- 0.58 # Based on Figure 8 in Reisen, W. K., Fang, Y., & Martinez, V. M. (2006). Effects of temperature on the transmission of west nile virus by Culex tarsalis (Diptera: Culicidae). Journal of medical entomology, 43(2), 309–17. Retrieved from http://www.ncbi.nlm.nih.gov/pubmed/16619616
  b_L1 <<- 1/(5*365) # birth rate for livestock
  mu_L1 <<- 1/(5*365) # death rate for livestock
  q_divided <<- 1E11 # parameter for the impact of the livestock on vector fecundity and gonotrophic cycles
  q_A <<- 0.007 # probability of transovarial transmission
  epsilon_L1 <<- 2/7# input by Sophie 1/2 - should it be 2/7?
  gamma_L1 <<- 1/30 # input by Sophie 1/3     - should it be 30?
  t_dep <<- 0.229 # average time for egg deposition = 0.229 days in lab ref(8) GL2018 supp
  
  days <- day_length
  
  # In order to make sure this file is consistent with master file:
  
  if (modeloutput == "allMozadultsmeanlast28days") {

    
    # Load packages:
    library(deSolve)
    library(rootSolve)
    library(tidyverse)
    library(gridExtra)
    
    # Functions required:
    source("Functions/system_of_RVF_equations_DEN.R", chdir = TRUE)
    source("Functions/development_rate.R", chdir = TRUE)
    source("Functions/mortality_rate.R", chdir = TRUE)
    source("Functions/oviposition_rate.R", chdir = TRUE)
    source("Functions/Sharpe_DeMichele_parameters.R", chdir=TRUE)
    
    # SET PARAMETERS FROM ARGUMENT--------------------------------------------------
    
    assign("lambda_A", value = lambda_A, pos = 2, inherits = TRUE)
    assign("rho_max_A", value = rho_max_A, pos = 2, inherits = TRUE)
    assign("kappa_A", value = kappa_A, pos = 2, inherits = TRUE)
    assign("A_A", value = A_A, pos = 2, inherits = TRUE)
    
    # Temperature, Water Body Area (wbarea), Livestock Density - Set up
    
    database <- data.frame(daynum = 1:5000, tmean = constanttemp, wbarea1 = constantwb,
                           cellarea = cellarea, livestockda = livestocktotal) # this is not important in this coded/utilized
    
    # propsearcharea represents the fraction of the total cell area that mosquitoes are likely to search in - they will search the whole cell
    database$propsearcharea <- 1
    # wbarea represents the portion of the water body area that mosquitoes are likely to interact with.
    database$wbarea <- database$wbarea1 * database$propsearcharea
    
    # RUN Aedes SIMULATION ---------------------------------------------------------------
    
    # Only look at relevant data for the location(s) being studied - this is taken from the original code, but in essence it is the same spot each time
    database_thisloc <- database
    
    assign("database_thisloc", value = database_thisloc, pos = 2, inherits = TRUE)
    
    # Sequence from 0, 1, 2, 3,... until last day number:
    times <- seq(0, days, by = delta_t)
    
    # Additional functions:
    source("Functions/ponds.R", chdir = TRUE) 
    source("Functions/biting_rate.R", chdir = TRUE)
    
    # This function is not needed for this code as the temperature is set as a constant - this has been addressed in the system of RVF equations file.
    # source("Functions/temperature.R", chdir = TRUE) 

    # State Variables - calculated using only Aedes with desiccation =0 and delta = 1:
    # State_RVF <- c(O_a1 = 4380798686,
    #                O_a2 = 0,
    #                L_a = 15601261343,
    #                P_a = 5676180447,
    #                A1_a = 2861290378,
    #                F_a = 1448833708,
    #                A2_a = 9638652449)
    
    # Used to obtain state variables:
    State_RVF <- c(O_a1 = 100,
                   O_a2 = 0,
                   L_a = 0,
                   P_a = 0,
                   A1_a = 0,
                   F_a = 0,
                   A2_a = 0)
    

    N_L1 <- livestocktotal
    
    # Parameters to feed into ode function:
    Pars_RVF <- c(q_A, 
                  rho_max_A,
                  b_L1, 
                  mu_L1, 
                  epsilon_L1, 
                  gamma_L1, 
                  q_divided, 
                  as.numeric(N_L1), 
                  prop_find,
                  livestocktotal,
                  lambda_A,
                  A_A,
                  kappa_A
    )
    
    # Function for differential equations:
    out <- ode(y = State_RVF,
               times = times,
               func = system_of_RVF_equations_density_temp_Aedes,
               parms = Pars_RVF,
               atol = 1e-1, rtol = 1e-1, maxsteps = 500000
               )
    
    out_RVF_limiting <<- as.data.frame(out)
    
    # Add results from this loc to vector:
    
    cell_ode_output <- out_RVF_limiting 
    
    # Take data from all locs and bind into single dataframe:
    data_output <- cell_ode_output
    
    colnames(data_output) <-
      c("Time",
        "Immature Eggs",               
        "Mature Eggs",                
        "Larvae",                     
        "Pupae",                      
        "Nulliparous",        
        "Flyers",               
        "Parous")      
    
    facetPlotData <- data_output %>%
      tidyr::pivot_longer(cols=c("Immature Eggs",               
                                 "Mature Eggs",
                                 "Larvae",                     
                                 "Pupae",                      
                                 "Nulliparous",          
                                 "Flyers",                     
                                 "Parous"),                       
                          names_to = "Compartment",
                          values_to = "Value")
    
    facetPlotData$Compartment <- fct_inorder(facetPlotData$Compartment) # to keep them in order as above for plotting
    
    adults_cols <- c("Nulliparous",          
                     "Flyers",                     
                     "Parous")       
    
    toReturn <<- data_output
    
    toReturn$allMozadults <-rowSums(toReturn[, adults_cols])
    
    # if chosen "one", then the return plot will be the state variable plot and show the population. 
    if (what == "one") {
      facetPlot <- ggplot(data = facetPlotData, aes(x = Time, y = Value)) +
        geom_line(linewidth = 0.8) +
        facet_wrap(~Compartment, scales = "free") +
        theme_classic(base_size = 12) +
        labs(x = "Time (days)", y = "Population") +
        scale_y_continuous(
          expand = expansion(mult = c(0.1, 0.05)),  # Add 10% padding to the lower limit
          breaks = function(limits) {
            # Generate pretty breaks and ensure the lower limit is slightly above 0
            pretty_breaks <- scales::pretty_breaks(n = 5)(limits)
            pretty_breaks[pretty_breaks >= max(min(limits), 1e-6)]  # Exclude 0 and negatives
          }
        ) +
        theme(
          # plot.title = element_text(hjust = 0.5), strip.text = element_text(size = 7)
          plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),  # Centered and prominent title
          axis.title = element_text(size = 10),      # Axis titles are clear and legible
          axis.text = element_text(size = 8),       # Tick labels are a bit smaller
          strip.text = element_text(size = 7, face = "bold"),  # Facet labels are bold for clarity
          panel.grid.major = element_line(color = "gray90"),   # Optional: add light grid lines for reference
          panel.grid.minor = element_blank(),        # Remove minor grid lines to avoid clutter
          panel.spacing = unit(0.5, "lines")           # Adjust spacing between facets
        )

      facetPlot2 <- ggplot(data = facetPlotData, aes(x = Time, y = Value)) +
        geom_line(linewidth = 0.8) +
        facet_grid(~Compartment, scales = "fixed") +
        theme_classic(base_size = 12) +
        labs(x = "Time (days)", y = "Population") +
        scale_y_continuous(
          breaks = seq(0, 2e11, by = 2.5e10),  # Manually specify evenly spaced breaks
          labels = scales::label_number(scale = 1, accuracy = 1e10)  # Format labels for clarity
        ) +
        theme(
          plot.title = element_text(hjust = 0.5, size = 8, face = "bold"),  # Centered and prominent title
          axis.title = element_text(size = 8),      # Axis titles are clear and legible
          axis.text = element_text(size = 6),       # Tick labels are a bit smaller
          strip.text = element_text(size = 5, face = "bold"),  # Facet labels are bold for clarity
          panel.grid.major = element_line(color = "gray90"),   # Optional: add light grid lines for reference
          panel.grid.minor = element_blank(),        # Remove minor grid lines to avoid clutter
          panel.spacing = unit(1, "lines")           # Adjust spacing between facets
        )
      facetPlot2

      ggsave(plot = facetPlot,
             filename = paste0("Output/", "Aedes_SingleSim",  "_facetPlot.png"),
             width = 6, height = 4)
      ggsave(plot = facetPlot,
             filename = paste0("Output/", "Aedes_SingleSim",  "_facetPlot.pdf"),
             width = 6, height = 3,device=pdf)

      ggsave(plot = facetPlot2,
             filename = paste0("Output/", "Aedes_SingleSim",  "_facetPlot2.png"),
             width = 6, height = 3)
      ggsave(plot = facetPlot2,
             filename = paste0("Output/", "Aedes_SingleSim",  "_facetPlot2.pdf"),
             width = 6, height = 2,device=pdf)

      write.csv(data_output, file = paste0("Output/", "Aedes_SingleSim", ".csv"))
      
    }
    
    # If chosen "SA" then the return is the sum of all adults.
    if (what == "SA") {
      
      startday <- day_length-365
      mozMeanLast4weeks <- toReturn %>%
        dplyr::filter(Time > startday)
      return(sum(mozMeanLast4weeks$allMozadults))
    }
    
  }
  
  else {stop("Error - Model output in SA_SingleSimulation.R does not match mastercode")}
}