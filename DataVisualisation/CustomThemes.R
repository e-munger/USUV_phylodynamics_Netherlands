# Custom ggplot2 themes used to give a consistent look across figures
# (colors, font sizes, axis/legend styling). Sourced by the main
# visualisation script. Some final styling adjustments were made
# afterwards in Illustrator and are not reflected here.

Arbo_theme_grid <- function() { 
  theme_light(base_family = "Helvetica", base_size = 10) +
    theme(
      plot.title = element_text(face = "plain", 
                                margin = margin(0, 0, 4, 0), color = "#252525"),
      plot.subtitle = element_text(size = 8.89, face = "plain", 
                                   margin = margin(0, 0, 3, 0), color = "#737373"),
      axis.title = element_text(size = 8.89, face = "plain", 
                                color = "#252525"),
      
      # Line weight for axes
      axis.ticks.x = element_line(color = "#231f20"),
      axis.ticks.y = element_blank(),
      
      axis.line.x = element_line(color = "#231f20"), # Adjust axis line weight and color
      axis.line.y = element_blank(), # No y axis line
      
      # Specify axis text colors
      axis.text.x = element_text(size = 8.89, color = "#231f20"), # X-axis text color
      axis.text.y = element_text(size = 8.89, color = "#231f20"), # Y-axis text color
      
      # Grid line settings
      panel.grid.major = element_line(color = "#d1d3d4"),  
      panel.grid.minor = element_line(color = "#d1d3d4"),
      
      # Remove plot frame
      panel.border = element_blank(),         # No border around the plot
      
      # Legend settings
      legend.title = element_text(size = 7.9, face = "plain", color = "#252525",
                                  margin = margin(0, 0, 2, 0)),  # Reduce margin to bring items closer
      legend.text = element_text(size = 7.9, color = "#737373"),
      
      # Adjust spacing between symbols and text in the legend
      legend.spacing = unit(0.05, "cm"),              # Controls spacing between legend items
      legend.key.size = unit(0.2, "cm"),             # Controls the size of the legend key (symbols)
      legend.margin = margin(3, 3, 3, 3)             # Controls overall margin of legend
    )
}

Arbo_theme_axes <- function() { 
  theme_light(base_family = "Helvetica", base_size = 10) +
    theme(
      plot.title = element_text(face = "plain", 
                                margin = margin(0, 0, 4, 0), color = "#252525"),
      plot.subtitle = element_text(size = 8.89, face = "plain", 
                                   margin = margin(0, 0, 3, 0), color = "#737373"),
      axis.title = element_text(size = 8.89, face = "plain", 
                                color = "#252525"),
      
      # Line weight for axes
      axis.ticks.x = element_line(color = "#231f20"),
      axis.ticks.y = element_line(color = "#231f20"),
      
      axis.line.x = element_line(color = "#231f20"), # Adjust axis line weight and color
      axis.line.y = element_line(color = "#231f20"), # No y axis line
      
      # Specify axis text colors
      axis.text.x = element_text(size = 8.89, color = "#231f20", angle = 45, hjust = 1), # X-axis text color
      axis.text.y = element_text(size = 8.89, color = "#231f20"), # Y-axis text color
      
      # Grid line settings
      panel.grid.major = element_blank(),     # No major grid lines
      panel.grid.minor = element_blank(),     # No minor grid lines
      
      # Remove plot frame
      panel.border = element_blank(),         # No border around the plot
      
      # Legend settings
      legend.title = element_text(size = 7.9, face = "plain", color = "#252525",
                                  margin = margin(0, 0, 2, 0)),  # Reduce margin to bring items closer
      legend.text = element_text(size = 7.9, color = "#737373"),
      
      # Adjust spacing between symbols and text in the legend
      legend.spacing = unit(0.05, "cm"),              # Controls spacing between legend items
      legend.key.size = unit(0.2, "cm"),             # Controls the size of the legend key (symbols)
      legend.margin = margin(3, 3, 3, 3)             # Controls overall margin of legend
    )
}
