#Clear Workspace 
rm(list=ls())

#Load packages 
library(ggplot2)
library(tidyverse)

#Load royle-nichols models from rn_model_files subfolder in data folder

bacs_mod <- readRDS("data/rn_model_files/rn_model_t_bacs.rds")
blgr_mod <- readRDS("data/rn_model_files/rn_model_t_blgr.rds")
coni_mod <- readRDS("data/rn_model_files/rn_model_t_coni.rds")
coye_mod <- readRDS("data/rn_model_files/rn_model_t_coye.rds")
cwwi_mod <- readRDS("data/rn_model_files/rn_model_t_cwwi.rds")
inbu_mod <- readRDS("data/rn_model_files/rn_model_t_inbu.rds")
nobo_mod <- readRDS("data/rn_model_files/rn_model_t_nobo.rds")
praw_mod <- readRDS("data/rn_model_files/rn_model_t_praw.rds")
wevi_mod <- readRDS("data/rn_model_files/rn_model_t_wevi.rds")

