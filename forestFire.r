library(magick)

# PARAMETERS
grid_size <- 50
tree_density <- 0.6
p_spread_base <- 0.3  # Increased for better spread
wind_direction <- "E"  # "N", "S", "E", "W", or "NONE"
wind_boost <- 0.3
heat_threshold <- 2   # Number of burning neighbors that guarantees ignition
n_steps <- 200
frame_dir <- "frames"
gif_path <- "forest_fire.gif"
png_dir <- "png_frames"  # Directory to save individual PNG frames

# STATE CODES
EMPTY <- 0
TREE <- 1
BURNING <- 2
BURNED <- 3  # New state to track burned cells separately from initially empty cells

dir.create(frame_dir, showWarnings = FALSE)
dir.create(png_dir, showWarnings = FALSE)

# INIT FUNCTIONS
initialize_forest <- function(size, density) {
  forest <- matrix(EMPTY, nrow = size, ncol = size)
  forest[] <- ifelse(runif(size * size) < density, TREE, EMPTY)
  return(forest)
}

ignite_random_tree <- function(forest) {
  tree_positions <- which(forest == TREE, arr.ind = TRUE)
  if (nrow(tree_positions) == 0) return(forest)
  
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
  return(forest)
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

# STEP FIRE function with enhanced neighborhood effects
step_fire <- function(forest, p_spread, wind_direction, wind_boost, heat_threshold) {
  size <- nrow(forest)
  new_forest <- forest
  fire_spread_cells <- which(forest == BURNING, arr.ind = TRUE)
  
  # Track the number of burning neighbors for each tree
  burning_neighbors <- matrix(0, nrow = size, ncol = size)
  
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
        burning_neighbors[ni, nj] <- burning_neighbors[ni, nj] + 1
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
          next
        }
        
        # Otherwise, use probability-based spread
        prob <- p_spread  # Base spread probability
        
        # Increase probability based on number of burning neighbors
        prob <- prob + (burning_neighbors[ni, nj] * 0.1)
        
        # Apply wind boost in the wind direction
        if (dir == wind_direction) {
          prob <- prob + wind_boost
        } else if (wind_direction == "N" && (dir == "NE" || dir == "NW")) {
          prob <- prob + wind_boost * 0.5
        } else if (wind_direction == "S" && (dir == "SE" || dir == "SW")) {
          prob <- prob + wind_boost * 0.5
        } else if (wind_direction == "E" && (dir == "NE" || dir == "SE")) {
          prob <- prob + wind_boost * 0.5
        } else if (wind_direction == "W" && (dir == "NW" || dir == "SW")) {
          prob <- prob + wind_boost * 0.5
        }
        
        # Cap probability at 0.95 to maintain some randomness
        prob <- min(prob, 0.95)
        
        # If the random chance is below the probability, ignite the tree
        if (runif(1) < prob) {
          new_forest[ni, nj] <- BURNING
        }
      }
    }
  }
  
  # Keep some trees burning for another step (30% chance)
  # Reduced from previous 50% for faster transition through burning phase
  for (i in 1:nrow(fire_spread_cells)) {
    x <- fire_spread_cells[i, 1]
    y <- fire_spread_cells[i, 2]
    
    # 70% chance to burn out, 30% chance to keep burning
    if (runif(1) > 0.3) {
      new_forest[x, y] <- BURNED
    }
    # else: keep it as BURNING
  }
  
  return(new_forest)
}

# PLOT AND SAVE FRAME
plot_and_save_forest <- function(forest, step, burned_percent) {
  # Save each frame as a PNG file
  png_filename <- sprintf("%s/frame_%03d.png", png_dir, step)
  png(png_filename, width = 500, height = 500)
  # Use a 4-color palette: gray for empty, green for trees, red for burning, dark gray for burned
  image(t(apply(forest, 2, rev)),
        col = c("gray", "forestgreen", "red", "darkgray"),
        axes = FALSE,
        main = paste("Step", step, "- Burned:", round(burned_percent, 1), "%"))
  dev.off()
  
  # Optionally, save the same frame to the frame directory (for GIF)
  frame_filename <- sprintf("%s/frame_%03d.png", frame_dir, step)
  png(frame_filename, width = 500, height = 500)
  image(t(apply(forest, 2, rev)),
        col = c("gray", "forestgreen", "red", "darkgray"),
        axes = FALSE,
        main = paste("Step", step, "- Burned:", round(burned_percent, 1), "%"))
  dev.off()
}

# RUN SIMULATION
forest <- initialize_forest(grid_size, tree_density)
initial_empty_count <- sum(forest == EMPTY)

# Ignite multiple trees for a more robust start
forest <- ignite_random_tree(forest)
# Ignite a few more trees (optional - more aggressive start)
for (i in 1:3) {
  forest <- ignite_random_tree(forest)
}

initial_tree_count <- sum(forest == TREE) + sum(forest == BURNING)
steps_taken <- 0

for (step in 1:n_steps) {
  steps_taken <- step
  
  # Calculate burned percentage based on both burned and burning trees
  burned_count <- sum(forest == BURNED)
  burning_count <- sum(forest == BURNING)
  burned_percent <- 100 * (burned_count + burning_count) / initial_tree_count
  
  plot_and_save_forest(forest, step, burned_percent)
  
  # Enhanced termination condition: either no burning trees or reached high burn percentage
  if (!any(forest == BURNING) || burned_percent > 99) break
  
  # Pass the heat_threshold parameter to the step_fire function
  forest <- step_fire(forest, p_spread_base, wind_direction, wind_boost, heat_threshold)
}

# CREATE GIF
frame_files <- list.files(frame_dir, pattern = "png$", full.names = TRUE)
animation <- image_read(frame_files)
animation <- image_animate(animation, fps = 5)  # Slow down for better visualization
image_write(animation, gif_path)

cat("GIF saved as:", gif_path, "\n")