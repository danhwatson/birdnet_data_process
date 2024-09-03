library(tidyverse)
library(glue)
library(lubridate)

# Data from all weather NOAA stations 
inventory_url <- "https://www.ncei.noaa.gov/pub/data/ghcn/daily/ghcnd-inventory.txt"

# Make a df for it 
inventory <- read_table(inventory_url,
                        col_names = c("station", "lat", "lon", "variable", "start", "end"))

# Read in to find site coordinates
site_cords <- read_csv("data/aru_site_coordinates.csv")

#creating values for lat and lon for each site treatment 
chemours_lat <- 31.03901 * 2 * pi / 360
chemours_lon <- -81.96080 * 2 * pi / 360

sansa_lat <-  31.44399 * 2 * pi / 360
sansa_lon <- --81.67899 * 2 * pi / 360

okefenokee_lat <- 30.69881 * 2 * pi / 360
okefenokee_lon <- -82.16504 * 2 * pi / 360

# Distance, d = 3963.0 * arccos[(sin(lat1) * sin(lat2)) + cos(lat1) * cos(lat2) * cos(long2 – long1)]
# The obtained distance, d, is in miles. If you want your value to be in units of kilometers, multiple d by 1.609344.
# d in kilometers = 1.609344 * d in miles
my_station_chemours <- inventory %>%
  mutate(lat_r = lat *2 *pi/360,
         lon_r = lon *2 *pi/360,
         d = 1.609344 * 3963 * acos((sin(lat_r) * sin(chemours_lat)) + cos(lat_r) * cos(chemours_lat) * cos(chemours_lon - lon_r))
  ) %>%
  filter(start == 2024) %>% #filter time frame
  top_n(n = -1, d) %>% #find the closest station
  distinct(station) %>% #remove duplicates
  pull(station)

my_station_sansa <- inventory %>% 
  mutate(lat_r = lat *2 *pi/360,
         lon_r = lon *2 *pi/360,
         d = 1.609344 * 3963 * acos((sin(lat_r) * sin(sansa_lat)) + cos(lat_r) * cos(sansa_lat) * cos(sansa_lon - lon_r))
  ) %>%
  filter(start == 2024) %>%
  top_n(n = -1, d) %>%
  distinct(station) %>%
  pull(station)

my_station_okefenokee <- inventory %>% 
  mutate(lat_r = lat *2 *pi/360,
         lon_r = lon *2 *pi/360,
         d = 1.609344 * 3963 * acos((sin(lat_r) * sin(okefenokee_lat)) + cos(lat_r) * cos(okefenokee_lat) * cos(okefenokee_lon - lon_r))
  ) %>%
  filter(start == 2024) %>%
  top_n(n = -1, d) %>%
  distinct(station) %>%
  pull(station)


# Read in data from each station
station_daily_chemours <- glue("https://www.ncei.noaa.gov/pub/data/ghcn/daily/by_station/USC00099502.csv.gz")
station_daily_sansa <- glue("https://www.ncei.noaa.gov/pub/data/ghcn/daily/by_station/LYM00062012.csv.gz")
station_daily_okefenokee <- glue("https://www.ncei.noaa.gov/pub/data/ghcn/daily/by_station/USC00099502.csv.gz")


local_weather_chemours <- read_csv(station_daily_chemours,
                          col_names = c("station", "date", "variable", "value", "a", "b", "c", "d")) %>%
  select(date, variable, value) %>%
  pivot_wider(names_from = "variable", values_from="value",
              values_fill = 0) %>%
  mutate(date = ymd(date))
