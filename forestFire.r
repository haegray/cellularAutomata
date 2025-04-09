library(magick)

# PARAMETERS
grid_size <- 50
tree_density <- 0.6
p_spread_base <- 0.3  # Base probability for fire spread
wind_direction <- sample(c("N", "S", "E", "W", "NE", "NW", "SE", "SW", "NONE"), 1)  # Random wind direction
wind_strength <- sample(c("none", "light", "moderate", "strong"), 1, 
                        prob = c(0.1, 0.3, 0.4, 0.2))  # Random wind strength

# Set wind boost based on strength
if (wind_strength == "none") {
  wind_boost <- 0
  wind_direction <- "NONE"
} else if (wind_strength == "light") {
  wind_boost <- 0.15
} else if (wind_strength == "moderate") {
  wind_boost <- 0.3
} else { # strong
  wind_boost <- 0.5
}

cat("Wind direction:", wind_direction, "- Wind strength:", wind_strength, 
    "(boost:", wind_boost, ")\n")

heat_threshold <- 2   # Number of burning neighbors that guarantees ignition
max_burn_time <- 3    # Maximum steps a tree can burn before being consumed
persistence_factor <- 0.7  # Chance of a tree continuing to burn based on neighbors
n_steps <- 200
frame_dir <- "frames"
gif_path <- "forest_fire.gif"
png_dir <- "png_frames"  # Directory to save individual PNG frames

# STATE CODES
EMPTY <- 0
TREE <- 1
BURNING <- 2  # Now will have associated burn_time value
BURNED <- 3  # Fully consumed by fire

# Create directories
dir.create(frame_dir, showWarnings = FALSE)
dir.create(png_dir, showWarnings = FALSE)

# INIT FUNCTIONS
initialize_forest <- function(size, density) {
  forest <- matrix(EMPTY, nrow = size, ncol = size)
  forest[] <- ifelse(runif(size * size) < density, TREE, EMPTY)
  return(forest)
}

# Function to initialize burn time tracking matrix
initialize_burn_times <- function(size) {
  burn_times <- matrix(0, nrow = size, ncol = size)
  return(burn_times)
}

ignite_random_tree <- function(forest, burn_times) {
  tree_positions <- which(forest == TREE, arr.ind = TRUE)
  if (nrow(tree_positions) == 0) return(list(forest = forest, burn_times = burn_times))
  
  # Select a random tree from the center 50% of the grid for better spread
  size <- nrow(forest)
  center_range <- floor(size * 0.25):(ceiling(size * 0.75))
  
  # Filter tree positions to favor those in the center area
  center_trees <- tree_positions[
    tree_positions[,1] %in% center_range & 
    tree_positions[,2] %in% center_range, 
  ]
  
  # If there are center trees, choose from them, otherwise choose from all trees
  if (nrow(center_trees) > 0) {
    chosen <- center_trees[sample(nrow(center_trees), 1), ]
  } else {
    chosen <- tree_positions[sample(nrow(tree_positions), 1), ]
  }
  
  forest[chosen[1], chosen[2]] <- BURNING
  burn_times[chosen[1], chosen[2]] <- 1  # Start burn timer at 1
  
  return(list(forest = forest, burn_times = burn_times))
}

# Define direction mappings
direction_map <- list(
  "N" = c(0, 1),   # North
  "S" = c(0, -1),  # South
  "E" = c(1, 0),   # East
  "W" = c(-1, 0),  # West
  "NE" = c(1, 1),  # Northeast
  "NW" = c(-1, 1), # Northwest
  "SE" = c(1, -1), # Southeast
  "SW" = c(-1, -1) # Southwest
)

# Simplified neighbor function that returns coordinates directly
get_neighbors <- function(i, j, size) {
  directions <- expand.grid(dx = -1:1, dy = -1:1)
  directions <- directions[!(directions$dx == 0 & directions$dy == 0), ]
  neighbors <- list()
  
  for (k in 1:nrow(directions)) {
    dx <- directions$dx[k]
    dy <- directions$dy[k]
    ni <- i + dx
    nj <- j + dy
    
    if (ni >= 1 && ni <= size && nj >= 1 && nj <= size) {
      # Store the neighbor coordinates along with its direction
      dir_key <- "NONE"
      if (dx == 1 && dy == 0) dir_key <- "E"
      else if (dx == -1 && dy == 0) dir_key <- "W"
      else if (dx == 0 && dy == 1) dir_key <- "N"
      else if (dx == 0 && dy == -1) dir_key <- "S"
      else if (dx == 1 && dy == 1) dir_key <- "NE"
      else if (dx == -1 && dy == 1) dir_key <- "NW"
      else if (dx == 1 && dy == -1) dir_key <- "SE"
      else if (dx == -1 && dy == -1) dir_key <- "SW"
      
      neighbors[[length(neighbors) + 1]] <- list(
        coords = c(ni, nj),
        direction = dir_key
      )
    }
  }
  return(neighbors)
}

# Enhanced STEP FIRE function with burn time and persistence
step_fire <- function(forest, burn_times, p_spread, wind_direction, wind_boost, heat_threshold, 
                      max_burn_time, persistence_factor) {
  size <- nrow(forest)
  new_forest <- forest
  new_burn_times <- burn_times
  
  fire_spread_cells <- which(forest == BURNING, arr.ind = TRUE)
  
  # If no burning cells, just return the current state
  if (nrow(fire_spread_cells) == 0) {
    return(list(forest = new_forest, burn_times = new_burn_times))
  }
  
  # Track the number of burning neighbors for each tree
  burning_neighbors <- matrix(0, nrow = size, ncol = size)
  
  # Calculate burning intensity for each burning tree (1 to max_burn_time)
  burning_intensity <- matrix(0, nrow = size, ncol = size)
  for (i in 1:nrow(fire_spread_cells)) {
    x <- fire_spread_cells[i, 1]
    y <- fire_spread_cells[i, 2]
    burning_intensity[x, y] <- burn_times[x, y] / max_burn_time
  }
  
  # Count burning neighbors for all cells
  for (i in 1:nrow(fire_spread_cells)) {
    x <- fire_spread_cells[i, 1]
    y <- fire_spread_cells[i, 2]
    
    neighbors <- get_neighbors(x, y, size)
    for (n in neighbors) {
      ni <- n$coords[1]
      nj <- n$coords[2]
      
      # Increment the burning neighbor count for this cell
      if (forest[ni, nj] == TREE) {
        # Weight burning neighbors by their intensity
        intensity_factor <- burn_times[x, y] / max_burn_time
        burning_neighbors[ni, nj] <- burning_neighbors[ni, nj] + intensity_factor
      }
    }
  }
  
  # Process fire spread with probability and heat threshold logic
  for (i in 1:nrow(fire_spread_cells)) {
    x <- fire_spread_cells[i, 1]
    y <- fire_spread_cells[i, 2]
    
    neighbors <- get_neighbors(x, y, size)
    for (n in neighbors) {
      ni <- n$coords[1]
      nj <- n$coords[2]
      dir <- n$direction
      
      if (forest[ni, nj] == TREE) {
        # If tree has enough burning neighbors, it will definitely catch fire
        if (burning_neighbors[ni, nj] >= heat_threshold) {
          new_forest[ni, nj] <- BURNING
          new_burn_times[ni, nj] <- 1  # Start burn timer
          next
        }
        
        # Otherwise, use probability-based spread
        prob <- p_spread  # Base spread probability
        
        # Increase probability based on number of burning neighbors
        prob <- prob + (burning_neighbors[ni, nj] * 0.1)
        
        # Apply wind boost for various wind directions
        # Primary directions
        if (dir == wind_direction) {
          prob <- prob + wind_boost
        }
        # Handle diagonal wind directions
        else if (wind_direction == "NE" && (dir == "N" || dir == "E")) {
          prob <- prob + wind_boost * 0.7
        }
        else if (wind_direction == "NW" && (dir == "N" || dir == "W")) {
          prob <- prob + wind_boost * 0.7
        }
        else if (wind_direction == "SE" && (dir == "S" || dir == "E")) {
          prob <- prob + wind_boost * 0.7
        }
        else if (wind_direction == "SW" && (dir == "S" || dir == "W")) {
          prob <- prob + wind_boost * 0.7
        }
        # Adjacent to wind direction (primary directions)
        else if (wind_direction == "N" && (dir == "NE" || dir == "NW")) {
          prob <- prob + wind_boost * 0.5
        }
        else if (wind_direction == "S" && (dir == "SE" || dir == "SW")) {
          prob <- prob + wind_boost * 0.5
        }
        else if (wind_direction == "E" && (dir == "NE" || dir == "SE")) {
          prob <- prob + wind_boost * 0.5
        }
        else if (wind_direction == "W" && (dir == "NW" || dir == "SW")) {
          prob <- prob + wind_boost * 0.5
        }
        
        # Cap probability at 0.95 to maintain some randomness
        prob <- min(prob, 0.95)
        
        # If the random chance is below the probability, ignite the tree
        if (runif(1) < prob) {
          new_forest[ni, nj] <- BURNING
          new_burn_times[ni, nj] <- 1  # Start burn timer
        }
      }
    }
  }
  
  # Process burning trees for persistence or burnout
  for (i in 1:nrow(fire_spread_cells)) {
    x <- fire_spread_cells[i, 1]
    y <- fire_spread_cells[i, 2]
    
    # Calculate persistence chance based on neighboring burning trees
    neighbors <- get_neighbors(x, y, size)
    burning_neighbor_count <- 0
    for (n in neighbors) {
      ni <- n$coords[1]
      nj <- n$coords[2]
      if (forest[ni, nj] == BURNING) {
        burning_neighbor_count <- burning_neighbor_count + 1
      }
    }
    
    # Persistence chance increases with more burning neighbors
    # Base: persistence_factor, +5% for each burning neighbor
    persistence_chance <- persistence_factor + (0.05 * burning_neighbor_count)
    persistence_chance <- min(persistence_chance, 0.9)  # Cap at 90%
    
    # Increment burn time for this cell
    current_burn_time <- burn_times[x, y] + 1
    
    # If burn time exceeds max or fails persistence check, mark as BURNED
    if (current_burn_time > max_burn_time || runif(1) > persistence_chance) {
      new_forest[x, y] <- BURNED
      new_burn_times[x, y] <- 0
    } else {
      # Continue burning with incremented burn time
      new_burn_times[x, y] <- current_burn_time
    }
  }
  
  return(list(forest = new_forest, burn_times = new_burn_times))
}

# Enhanced PLOT AND SAVE function with burn intensity visualization
plot_and_save_forest <- function(forest, burn_times, step, burned_percent, max_burn_time) {
  # Create visualization matrix (for coloring)
  viz_matrix <- matrix(0, nrow = nrow(forest), ncol = ncol(forest))
  
  # Map states to visualization indices
  # 1 = Empty (gray)
  # 2 = Tree (green)
  # 3 = Burned (dark gray)
  # 4+ = Burning with different intensities (orange to dark red)
  
  # First set base states
  viz_matrix[forest == EMPTY] <- 1  # Empty
  viz_matrix[forest == TREE] <- 2   # Tree
  viz_matrix[forest == BURNED] <- 3 # Burned
  
  # Then set burning cells with their intensity levels
  burning_cells <- which(forest == BURNING, arr.ind = TRUE)
  if (nrow(burning_cells) > 0) {
    for (i in 1:nrow(burning_cells)) {
      x <- burning_cells[i, 1]
      y <- burning_cells[i, 2]
      burn_intensity <- burn_times[x, y]
      # Map intensity to color index (add offset to avoid conflict with other states)
      viz_matrix[x, y] <- 3 + burn_intensity
    }
  }
  
  # Create custom color palette for different burning intensities
  burning_colors <- colorRampPalette(c("orange", "red", "darkred"))(max_burn_time)
  
  # Create complete color palette
  col_palette <- c("gray", "forestgreen", "darkgray")
  col_palette <- c(col_palette, burning_colors)
  
  # Save each frame as a PNG file
  png_filename <- sprintf("%s/frame_%03d.png", png_dir, step)
  png(png_filename, width = 500, height = 500)
  
  image(t(apply(viz_matrix, 2, rev)),
        col = col_palette,
        axes = FALSE,
        main = paste("Step", step, "- Burned:", round(burned_percent, 1), "%"))
  
  # Add legend for burn intensity
  legend("topright", 
         legend = c("Empty", "Tree", "Burned", paste("Burning", 1:max_burn_time)), 
         fill = c("gray", "forestgreen", "darkgray", burning_colors),
         cex = 0.7, bg = "white")
  
  dev.off()
  
  # Optionally, save the same frame to the frame directory (for GIF)
  frame_filename <- sprintf("%s/frame_%03d.png", frame_dir, step)
  png(frame_filename, width = 500, height = 500)
  
  image(t(apply(viz_matrix, 2, rev)),
        col = col_palette,
        axes = FALSE,
        main = paste("Step", step, "- Burned:", round(burned_percent, 1), "%"))
  
  # Add legend for burn intensity
  legend("topright", 
         legend = c("Empty", "Tree", "Burned", paste("Burning", 1:max_burn_time)), 
         fill = c("gray", "forestgreen", "darkgray", burning_colors),
         cex = 0.7, bg = "white")
  
  dev.off()
}

# RUN SIMULATION
forest <- initialize_forest(grid_size, tree_density)
burn_times <- initialize_burn_times(grid_size)
initial_empty_count <- sum(forest == EMPTY)

# Ignite initial trees
result <- ignite_random_tree(forest, burn_times)
forest <- result$forest
burn_times <- result$burn_times

# Ignite a few more trees
for (i in 1:3) {
  result <- ignite_random_tree(forest, burn_times)
  forest <- result$forest
  burn_times <- result$burn_times
}

initial_tree_count <- sum(forest == TREE) + sum(forest == BURNING)
steps_taken <- 0

for (step in 1:n_steps) {
  steps_taken <- step
  
  # Calculate burned percentage based on both burned and burning trees
  burned_count <- sum(forest == BURNED)
  burning_count <- sum(forest == BURNING)
  burned_percent <- 100 * (burned_count + burning_count) / initial_tree_count
  
  plot_and_save_forest(forest, burn_times, step, burned_percent, max_burn_time)
  
  # Enhanced termination condition: either no burning trees or reached high burn percentage
  if (!any(forest == BURNING) || burned_percent > 99) break
  
  # Step the simulation
  result <- step_fire(forest, burn_times, p_spread_base, wind_direction, wind_boost, 
                      heat_threshold, max_burn_time, persistence_factor)
  forest <- result$forest
  burn_times <- result$burn_times
}

# CREATE GIF
frame_files <- list.files(frame_dir, pattern = "png$", full.names = TRUE)
animation <- image_read(frame_files)
animation <- image_animate(animation, fps = 5)  # Slow down for better visualization
image_write(animation, gif_path)

cat("GIF saved as:", gif_path, "\n")