library(httr)
library(jsonlite)
library(dplyr)
library(ggplot2)
setwd("C:/Users/pm83056/OneDrive - Office National des Forets/Bureau/Spatial_Buxbaumia_viridis/01_Statistics")

# Fetch data from OpenAlex API
url <- "https://api.openalex.org/works?filter=title.search:Buxbaumia%20viridis&per-page=200"
res <- fromJSON(content(GET(url), "text"), flatten = TRUE)

# Extract and clean publication years
df <- res$results %>%
  select(publication_year) %>%
  filter(!is.na(publication_year))

# Count publications per year
df_year <- df %>%
  group_by(publication_year) %>%
  summarise(n = n())

# Plot
ggplot(df_year, aes(x = publication_year, y = n)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(
    x = "Year",
    y = "Number of publications",
    title = "Publications containing 'Buxbaumia viridis' in the title, per year"
  )

# Save plot
#output_dir <- "02_Displays/Figures"

#ggsave(
#  filename = file.path(output_dir, "B_viridis_publications.png"),
#  width = 8, height = 5, dpi = 300)


