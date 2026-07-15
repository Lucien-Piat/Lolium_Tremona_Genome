library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthhires)
library(dplyr)
library(patchwork)
library(ggspatial)

df_coords <- read.csv("scripts/pop_gen/ch_coords.csv")
df_coords$lon <- ifelse(df_coords$lon > 100, -df_coords$lon, df_coords$lon)

q_matrix <- read.table("results/pca/Q_K6.tbl", skip = 1, header = FALSE)
colnames(q_matrix) <- c("accession", paste0("K", 1:6))

pop_map <- read.table("scripts/pop_gen/pop.tsv", header = FALSE, col.names = c("accession", "pop"))

q_pop <- merge(q_matrix, pop_map, by = "accession", all.x = TRUE)

q_pop$pop <- ifelse(grepl("^CH", q_pop$accession), q_pop$accession, q_pop$pop)

q_avg <- q_pop %>%
  group_by(pop) %>%
  summarise(across(starts_with("K"), mean)) %>%
  rename(accession = pop) 

df <- merge(df_coords, q_avg, by = "accession")

df$pop_group <- ifelse(grepl("^CH", df$accession), "CH", df$accession)

palette_df <- read.table("scripts/pop_gen/palete_ancestral.tsv", header = FALSE, comment.char = "")
cluster_colors <- palette_df$V1 
names(cluster_colors) <- paste0("K", 1:6)

point_palette_df <- read.table("scripts/pop_gen/palette_pop.tsv", header = FALSE, comment.char = "")
point_colors <- point_palette_df$V1
names(point_colors) <- point_palette_df$V2

df$lon_pie <- df$lon
df$lat_pie <- df$lat

# This is dirty but i gave up with this language 
df$lon_pie[df$accession == "GULF"] <- -124.2
df$lat_pie[df$accession == "GULF"] <- 44.0

df$lon_pie[df$accession == "L46"] <- -124.5
df$lat_pie[df$accession == "L46"] <- 45.4

df$lon_pie[df$accession == "L60"] <- -121.8
df$lat_pie[df$accession == "L60"] <- 44.3

df$lon_pie[df$accession == "L31"] <- -121.5
df$lat_pie[df$accession == "L31"] <- 45.7

df$lon_pie[df$accession == "PR"] <- -123.5
df$lat_pie[df$accession == "PR"] <- 40.2

df$lon_pie[df$accession == "SLB"] <- -120.5
df$lat_pie[df$accession == "SLB"] <- 38.8

df$lon_pie[df$accession == "CH4"] <- df$lon[df$accession == "CH4"] - 0.35
df$lat_pie[df$accession == "CH4"] <- df$lat[df$accession == "CH4"] + 0.25

df$lon_pie[df$accession == "TREM"] <- df$lon[df$accession == "TREM"] + 0.35
df$lat_pie[df$accession == "TREM"] <- df$lat[df$accession == "TREM"] - 0.20


build_map <- function(data_subset, map_bg, outer_border = NULL, pie_rad, crop_box = NULL) {
  p <- ggplot() +
    geom_sf(data = map_bg, fill = "white", color = "gray80", linewidth = 0.3)
  if (!is.null(outer_border)) {
    p <- p + geom_sf(data = outer_border, fill = NA, color = "black", linewidth = 0.9)
  }
  p <- p + geom_segment(data = data_subset, 
                        aes(x = lon, y = lat, xend = lon_pie, yend = lat_pie), 
                        color = "gray50", linewidth = 0.5)
  
  p <- p + geom_point(data = data_subset, aes(x = lon, y = lat, color = pop_group), size = 2.5) +
           scale_color_manual(values = point_colors)
  
  # Pie Charts
  for (i in 1:nrow(data_subset)) {
    pie_data <- data.frame(
      Cluster = paste0("K", 1:6),
      Proportion = as.numeric(data_subset[i, paste0("K", 1:6)])
    )
    
    pie_sub <- ggplot(pie_data, aes(x = "", y = Proportion, fill = Cluster)) +
      geom_col(width = 1) + 
      coord_polar("y", start = 0) +
      scale_fill_manual(values = cluster_colors) +
      theme_void() + theme(legend.position = "none")
    
    p <- p + annotation_custom(
      grob = ggplotGrob(pie_sub),
      xmin = data_subset$lon_pie[i] - pie_rad, xmax = data_subset$lon_pie[i] + pie_rad,
      ymin = data_subset$lat_pie[i] - pie_rad, ymax = data_subset$lat_pie[i] + pie_rad
    )
  }
  
  p <- p +
    geom_text(data = data_subset, aes(x = lon_pie, y = lat_pie, label = accession),
              vjust = -2.5, size = 3.5, fontface = "bold")
  
  p <- p + 
    annotation_scale(location = "bl", width_hint = 0.25, text_cex = 0.6) +
    theme_void() + 
    theme(
      legend.position = "none",
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1) 
    )
    
  if (!is.null(crop_box)) {
    p <- p + coord_sf(xlim = crop_box[1:2], ylim = crop_box[3:4], expand = FALSE)
  } else {
    p <- p + coord_sf(expand = FALSE)
  }
  
  return(p)
}


df_us <- df %>% filter(lon < 0)
df_ch <- df %>% filter(lon > 0)

us_states <- ne_states(country = "United States of America", returnclass = "sf") %>% 
  filter(name %in% c("California", "Oregon"))
us_border <- st_union(us_states)

swiss_cantons <- ne_states(country = "Switzerland", returnclass = "sf")
swiss_border <- ne_countries(country = "Switzerland", scale = "large", returnclass = "sf")

# Build Maps
map_us <- build_map(df_us, us_states, us_border, pie_rad = 0.35, crop_box = c(-125, -119, 38, 46.5))
map_ch <- build_map(df_ch, swiss_cantons, swiss_border, pie_rad = 0.12)

final_plot <- map_us + map_ch + 
  plot_layout(widths = c(0.75, 1)) & 
  theme(plot.margin = margin(0, 0, 0, 0, "cm"))

if (!dir.exists("results/map")) { dir.create("results/map", recursive = TRUE) }
ggsave("results/map/combined_map.pdf", plot = final_plot, width = 12, height = 7)