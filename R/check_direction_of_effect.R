



### CHECK DIRECTION OF EFFECT FOR REPLICATION OUTCOMES

library(dplyr)
library(readr)
library(stringr)


OUTCOME <- Sys.getenv("OUTCOME")
DISCOVERY_OUTCOME <- Sys.getenv("DISCOVERY_OUTCOME")


# read in discovery data
significant_res_discovery <- read_csv(str_c("full_results/significant_results_",DISCOVERY_OUTCOME,".txt"), col_types = cols())
significant_res_discovery1 <- significant_res_discovery[which((significant_res_discovery$method == "IVW" | significant_res_discovery$method == "Inverse variance weighted" | significant_res_discovery$method == "Wald ratio") & significant_res_discovery$clump_thresh == "0.2"),]
discovery <- significant_res_discovery1[,c("exposure","tissue","beta")]
names(discovery) <- c("exposure","tissue","discovery_beta")


# read in replication data
significant_res_replication <- read_csv(str_c("full_results/significant_results_",OUTCOME,".txt"), col_types = cols())
significant_res_replication1 <- significant_res_replication[which((significant_res_replication$method == "IVW" | significant_res_replication$method == "Inverse variance weighted" | significant_res_replication$method == "Wald ratio") & significant_res_replication$clump_thresh == "0.2"),]
replication <- significant_res_replication1[,c("exposure","tissue","beta")]
names(replication) <- c("exposure","tissue","replication_beta")


# put them into one data frame
direction <- right_join(discovery, replication, by = c("exposure","tissue"))
direction <- distinct(direction)
direction <- as.data.frame(direction)


# check if the direction of effect is consistent
direction$discovery_direction[direction$discovery_beta < 0]  <- "neg"
direction$discovery_direction[direction$discovery_beta > 0]  <- "pos"
direction$replication_direction[direction$replication_beta < 0]  <- "neg"
direction$replication_direction[direction$replication_beta > 0]  <- "pos"
direction$consistency <- "NO"
direction$consistency[direction$discovery_direction == direction$replication_direction] <- "YES"

# report to log file
if (all(direction$consistency == "YES")) {
  print("direction of effect is the SAME between discovery and replication for all outcomes")
} else if (all(direction$consistency == "NO")) {
  print("direction of effect is the OPPOSITE between discovery and replication for all outcomes")
} else {
  print("direction of effect is not the same for all outcomes")
}

# write out the results
write.table(direction, str_c("full_results/metrics_direction_of_effect_",OUTCOME,"_vs_",DISCOVERY_OUTCOME,".txt"))

print("mission complete")
