library("dplyr")
library("stringr")

# read in the previous file
trouble_already <- read.table("exposures_to_remove.txt", sep = ",", header = T, colClasses="character")
trouble_new <- data.frame()


# genes causing errors during clumping

trouble_new1 <- tryCatch({
        if (!is.na(file.size("clumping_trouble_exposures5.txt"))){
        trouble_clump <- read.table("clumping_trouble_exposures5.txt", header = F, colClasses="character")
        names(trouble_clump)[1] <- "exposure_data"
        names(trouble_clump)[4] <- "exposure"
        names(trouble_clump)[2] <- "outcome"
        trouble_new <- rbind(trouble_already,trouble_clump[,c("exposure_data","exposure","outcome")])
           }
        }, error = function(err) {
            # error handler picks up where error was generated
            print(paste("No more genes failed during clumping!:  ",err))
            return(data.frame())

        })



# genes causing errors during ld matrix method

trouble_new2 <- tryCatch({
        if (!is.na(file.size("ld_matrix_trouble_exposures2.txt"))){
        trouble_matrix <- read.table("ld_matrix_trouble_exposures2.txt", header = F, colClasses="character")
        names(trouble_matrix)[1] <- "exposure_data"
        names(trouble_matrix)[4] <- "exposure"
        names(trouble_matrix)[2] <- "outcome"
        trouble_new <- rbind(trouble_already,trouble_matrix[,c("exposure_data","exposure","outcome")])
           }
        }, error = function(err) {
            # error handler picks up where error was generated
            print(paste("No more genes failed during ld matrix method!:  ",err))
            return(data.frame())
        })

trouble_new <- rbind(trouble_new1, trouble_new2)

# write out if caught anything
if (plyr::empty(trouble_new) == TRUE) {
print("no more genes failed during any method")
} else {
write.table(distinct(trouble_new), "exposures_to_remove.txt", sep = ",", row.names = F)
}
