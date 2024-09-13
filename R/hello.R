#' Calibrate the GENEActiv Device
#'
#' This function reads the .bin file from a GENEActiv device and returns
#' the calibration scale and offset values.
#'
#' @param binfile The file path to the .bin file from the GENEActiv device.
#' @return A numeric vector containing the calibration scale and calibration offset.
#' @export
#' @examples
#' calib <- calibration("path/to/binfile.bin")
calibration <- function(binfile) {
  calibration_sphere <- GENEActiv.calibrate(binfile, printsummary = FALSE)
  calibration_offset <- calibration_sphere$offset
  calibration_scale <- calibration_sphere$scale
  return(c(calibration_scale, calibration_offset))
}

#' Get Data for a Specified Time Interval
#'
#' This function reads accelerometer data from the GENEActiv device for a specified time interval,
#' applies the calibration, and calculates the number of non-wear hours.
#'
#' @param v_bin_data The binary data read from the GENEActiv device.
#' @param binfile The file path to the .bin file from the GENEActiv device.
#' @param start The start time for the data extraction.
#' @param end The end time for the data extraction.
#' @param calib A calibration vector obtained from the `calibration()` function.
#' @param freq The frequency of measurements in Hz.
#' @return A list containing the following elements:
#' \item{DataFrame}{The extracted accelerometer data.}
#' \item{Date}{The date of the midpoint of the extracted data.}
#' \item{WeekDay}{The weekday of the extracted date.}
#' \item{DayNo}{The day number within the study.}
#' \item{NonWear}{The number of non-wear hours.}
#' @export
#' @examples
#' interval_data <- get_interval(v_bin_data, "path/to/binfile.bin", "2023-01-01", "2023-01-02", calib, 100)
get_interval <- function(v_bin_data, binfile, start, end, calib, freq) {
  day_data <- read.bin(binfile, start = start, end = end)$data.out
  day_data <- cbind(cbind(day_data[,1], ((calib[1:3] * day_data[,c('x', 'y', 'z')]) + calib[4:6])), day_data[,7])
  n_hours <- (nrow(day_data) / freq / 3600)
  non_wear_hours <- mean(rollapply(day_data, (10 * 60 * freq), non_wear_window, by = (freq * 60), by.column = FALSE)) * n_hours
  date <- as.Date(as.POSIXct(day_data[nrow(day_data) / 2, 1], origin = "1970-01-01"))
  day_of_week <- weekdays(date)
  return(list("DataFrame" = day_data, "Date" = date, "WeekDay" = day_of_week, "DayNo" = strsplit(start, " ")[[1]][1], "NonWear" = non_wear_hours))
}

#' Detect Non-Wear Windows
#'
#' This function identifies periods where the GENEActiv device was likely not worn based on
#' temperature and movement data.
#'
#' @param window A data frame containing accelerometer data and temperature values.
#' @return An integer status indicating non-wear (0) or wear (1).
#' @export
#' @examples
#' status <- non_wear_window(window_data)
non_wear_window <- function(window) {
  window2_start <- 0
  window2_end <- nrow(window) / 2
  window1_start <- window2_end + 1
  window1_end <- nrow(window)
  window1 <- window[window1_start:window1_end,]
  window2 <- window[window2_start:window2_end,]
  temp_in_window1 <- mean(window1[,5])
  temp_in_window2 <- mean(window2[,5])
  stds_in_window1 <- apply(window1[,c('x', 'y', 'z')], 2, sd)
  min_in_window1 <- apply(window1[,c('x', 'y', 'z')], 2, min)
  max_in_window1 <- apply(window1[,c('x', 'y', 'z')], 2, max)
  range_in_window1 <- max_in_window1 - min_in_window1
  t0 <- 26
  status <- 0
  stationairy <- range_in_window1 < 0.05 & stds_in_window1 < 0.013
  if ((sum(stationairy) >= 2) & (temp_in_window1 < t0)) {
    status <- 0
  } else if (temp_in_window1 > t0) {
    status <- 1
  } else {
    if (temp_in_window1 > temp_in_window2) {
      status <- 1
    }
    if (temp_in_window1 < temp_in_window2) {
      status <- 0
    }
  }
  return(status)
}

#' Add ENMO (Euclidean Norm Minus One) Metric
#'
#' This function computes the ENMO metric, a measure of physical activity intensity,
#' and adds it as a new column to the dataset.
#'
#' @param day_matrix A data frame containing accelerometer data, including 'x', 'y', and 'z' columns.
#' @return A data frame with the ENMO values added as a new column.
#' @export
#' @examples
#' day_matrix_with_enmo <- add_ENMO(day_matrix)
add_ENMO <- function(day_matrix) {
  ENMO <- pmax((((rowSums(day_matrix[,c('x', 'y', 'z')]^2))^(1/2)) - 1), 0)
  day_matrix <- cbind(day_matrix, ENMO)
  return(day_matrix)
}

#' Posture Detection Using SedUp
#'
#' This function detects periods of sedentary behavior or standing using a logistic regression model.
#'
#' @param y A numeric vector of accelerometer data (single axis).
#' @param fs The sampling frequency of the data.
#' @param option Integer specifying which posture detection option to use (1 or 2).
#' @return A numeric vector where 1 indicates standing and 0 indicates sedentary behavior.
#' @export
#' @examples
#' standing_periods <- SedUp(y_data, 100, 1)
SedUp <- function(y, fs, option) {
  y <- y * -1
  constant <- 10
  if (option == 1) {
    coeffs <- c(-2.390, 2.542, 42.394)
    uni_thr <- 0.362
    nSec <- 4
  } else if (option == 2) {
    coeffs <- c(-2.420, 2.616, 44.083)
    uni_thr <- 0.372
    nSec <- 6
  }
  y <- y[1:(floor(length(y) / fs / constant) * fs * constant)]
  sdv <- sapply(split(y, rep(1:(length(y) / fs), each = fs)), sd)
  window1 <- nSec * constant * fs / 2
  window2 <- round(nSec * constant / 2)
  medacc <- integer(length(sdv) / constant)
  medsd <- integer(length(sdv) / constant)
  j1 <- 0
  j2 <- 0
  for (i in seq(from = 1, to = length(sdv), by = constant)) {
    j1 <- j1 + 1
    i1_1 <- (i - 1) * fs + 1 - window1
    i1_2 <- (i - 1) * fs + 1 + constant * fs + window1
    b <- i1_1:i1_2
    medacc[j1] <- median(y[b[b > 0 & b < length(y)]], na.rm = TRUE)
    j2 <- j2 + 1
    i2_1 <- (i - 1) - window2
    i2_2 <- (i - 1) + 1 + constant + window2
    b <- i2_1:i2_2
    medsd[j2] <- median(sdv[b[b > 0 & b < length(sdv)]], na.rm = TRUE)
  }
  logit <- coeffs[1] + coeffs[2] * medacc + coeffs[3] * medsd
  phat <- exp(logit) / (1 + exp(logit))
  standing <- matrix(data = NA, nrow = length(phat), ncol = 1)
  standing[phat >= uni_thr] <- 1
  standing[phat < uni_thr] <- 0
  standing <- rep(standing, each = constant * fs)
  return(standing)
}

#' Get Raw Data for Bouts
#'
#' This function extracts raw accelerometer data for specified activity bouts.
#'
#' @param D A data frame containing accelerometer data.
#' @param bouts A data frame with start and end times for each bout.
#' @param freq The sampling frequency.
#' @return A list of data frames for each bout.
#' @export
get_raw_data <- function(D, bouts, freq) {
  return(apply(bouts, 1, function(x) { D[x[1]:x[2], ] }))
}

#' Extract Metrics from a Raw Bout
#'
#' This function extracts various metrics from a raw bout of activity, such as volume and intensity.
#'
#' @param raw_bout A data frame containing accelerometer data for a bout.
#' @param freq The sampling frequency.
#' @return A list of metrics for the bout.
#' @export
extract_metrics <- function(raw_bout, freq) {
  volume <- sum(raw_bout[,"ENMO"]) / freq
  intensity <- mean(raw_bout[,"ENMO"])
  duration <- nrow(raw_bout) / freq
  percentiles <- quantile(raw_bout[,'ENMO'], probs = seq(0, 1, 0.05))
  sd <- sd(raw_bout[,"ENMO"])
  wind_sd <- wind_SD(raw_bout[,"ENMO"])

  return(list(
    Volume = volume,
    Intensity = intensity,
    Duration = duration,
    Percentiles = percentiles,
    SD = sd,
    Windsorized_SD = wind_sd
  ))
}


#' Compute Bouts of Activity
#'
#' This function detects and computes periods of activity and inactivity based on accelerometer data.
#'
#' @param D A data frame containing accelerometer data.
#' @param freq The sampling frequency.
#' @param timestamps A vector of timestamps corresponding to the data.
#' @return A data frame with start and end times of each activity bout.
#' @export
#' @examples
#' bouts <- compute_bouts(D, 100, timestamps)
compute_bouts <- function(D, freq, timestamps) {
  X <- SedUp(D[,3], freq, 1)
  times <- as.numeric(format(as.POSIXct(timestamps, origin = "1970-01-01"), "%H"))
  X[times >= 0 & times < 6] <- 0
  rle_bouts <- rle(X)
  start <- head(cumsum(c(0, rle_bouts$lengths)), -1)
  end <- cumsum(rle_bouts$lengths)
  return(data.frame(start = start, end = end, lab = rle_bouts$values))
}

#' Compute Data for One Participant
#'
#' This function processes a single .bin file for a participant, computes posture changes, and saves results.
#'
#' @param file_name The name of the .bin file for the participant.
#' @param binfolder The path to the folder containing the .bin files.
#' @param outfolder The path to the output folder where results will be saved.
#' @return A list containing data frames for each day for that participant.
#' @export
#' @examples
#' person_dataframes <- compute_one_person("participant.bin", "path/to/binfolder", "path/to/outfolder")
compute_one_person <- function(file_name, binfolder, outfolder) {
  binfile <- file.path(binfolder, file_name)
  file_name_without_bin <- strsplit(file_name, '.bin')[[1]]
  specific_out_folder <- file.path(outfolder, file_name_without_bin)
  dir.create(specific_out_folder, recursive = TRUE)
  v_bin_data <- read.bin(binfile, virtual = TRUE)
  freq <- as.numeric(as.character(v_bin_data$header['Measurement_Frequency', 1]))
  calib <- calibration(binfile)
  time_stamp_list <- v_bin_data[["page.timestamps"]]
  d1 <- time_stamp_list[1]
  d2 <- time_stamp_list[length(time_stamp_list)]
  dates <- c(seq(d1, d2, by = "day"), d2)
  days <- paste(seq(length(dates)), "00:00:00")
  filenames <- paste(paste('Day', seq(length(dates)), sep = ''), ".csv", sep = "")

  # Create a list to store the data frames for each person
  person_dataframes <- list()

  for (i in (1:(length(filenames) - 1))) {
    start <- days[i]
    end <- days[i + 1]
    D1 <- get_interval(v_bin_data, binfile, start = start, end = end, calib = calib, freq = freq)
    if (nrow(D1$DataFrame) > 3000) {
      bouts <- compute_bouts(D1$DataFrame, freq, timestamps = D1$DataFrame[,1])
      ENMO_tab <- add_ENMO(D1$DataFrame)
      raw_bouts <- get_raw_data(ENMO_tab, bouts, freq)
      if (nrow(bouts) > 1) {
        list_of_bouts <- lapply(raw_bouts, extract_metrics, freq = freq)
        table_of_bouts <- t(do.call("cbind", list_of_bouts))
        rownames(table_of_bouts) <- NULL
        table_of_bouts <- cbind(table_of_bouts, "Date" = as.character(D1$Date))
        table_of_bouts <- cbind(table_of_bouts, "WeekDay" = as.character(D1$WeekDay))
        table_of_bouts <- cbind(table_of_bouts, "DayNo" = D1$DayNo)
        table_of_bouts <- cbind(table_of_bouts, "NonWear" = D1$NonWear)
        table_of_bouts <- cbind(table_of_bouts, "FILE" = file_name)
        table_of_bouts <- cbind(table_of_bouts, "Active" = bouts[,'lab'])
      }
      table_of_bouts <- na.omit(table_of_bouts)
      filepath <- file.path(specific_out_folder, filenames[i])
      write.csv(table_of_bouts, filepath, row.names = FALSE)
      # Store each day's dataframe in the list
      person_dataframes[[paste0("Day", i)]] <- table_of_bouts
    }
  }

  # Combine all days into a single file
  temp <- list.files(specific_out_folder, pattern = "*.csv")
  myfiles <- lapply(file.path(specific_out_folder, temp), read.csv)
  complete_file <- do.call('rbind', myfiles)
  write.csv(complete_file, file.path(specific_out_folder, "all_days.csv"), row.names = FALSE)

  # Return the list of dataframes for this person
  return(person_dataframes)
}


#' Compute Data for All Participants
#'
#' This function processes all .bin files in the specified folder, computes posture changes, and burstiness metrics.
#'
#' @param binfolder The path to the folder containing the .bin files.
#' @param outfolder The path to the output folder where results will be saved.
#' @return A list containing data frames for all participants.
#' @export
#' @examples
#' results <- compute_for_all("path/to/binfolder", "path/to/outfolder")
compute_for_all <- function(binfolder, outfolder) {
  all_binfiles <- list.files(binfolder, pattern = "*.bin")
  all_dataframes <- list()  # Store all dataframes from all participants

  for (i in all_binfiles) {
    person_dataframes <- compute_one_person(i, binfolder, outfolder)
    all_dataframes[[i]] <- person_dataframes  # Store each person's dataframes in the list
  }

  # Combine all individual 'all_days.csv' files into a single dataframe
  temp <- list.files(outfolder, pattern = "*all_days.csv", recursive = TRUE)
  myfiles <- lapply(file.path(outfolder, temp), read.csv)
  complete_file <- do.call('rbind', myfiles)
  write.csv(complete_file, file.path(outfolder, "all_people.csv"), row.names = FALSE)

  return(list("dataframes" = all_dataframes, "output_folder" = outfolder))
}


#' Compute Windsorized Standard Deviation
#'
#' This function calculates the Windsorized standard deviation, reducing the influence of outliers.
#'
#' @param vect A numeric vector of data.
#' @return The Windsorized standard deviation of the input vector.
#' @export
#' @examples
#' wind_sd <- wind_SD(data_vector)
wind_SD <- function(vect) {
  top <- as.numeric(quantile(vect, probs = c(0.95)))
  bottom <- as.numeric(quantile(vect, probs = c(0.05)))
  wind_vect <- pmin(pmax(vect, bottom), top)
  return(sd(wind_vect))
}

#' Calculate Novel Burstiness Metrics
#'
#' This function calculates the novel burstiness measure (An) and the burstiness measure
#' considering minimum interevent time (An_y) based on event timings.
#'
#' @param df A data frame containing 'Date', 'Time', and 'Active' columns.
#' @return A list of burstiness metrics, including novel burstiness, burstiness considering minimum interevent time, and related statistics.
#' @export
#' @examples
#' burstiness_metrics <- calculate_burstiness(df)
calculate_burstiness <- function(df) {
  df <- as.data.frame(df)
  required_columns <- c("Date", "Time", "Active")
  if (!all(required_columns %in% names(df))) {
    message("Skipping dataframe: Missing required columns.")
    return(NULL)
  }
  filtered_data <- df %>% filter(Active == 1)
  filtered_data$DateTime <- as.POSIXct(paste(filtered_data$Date, filtered_data$Time), format = "%Y-%m-%d %H:%M:%S")
  filtered_data <- filtered_data %>% distinct() %>%
    arrange(DateTime) %>%
    mutate(Time_Diff = as.numeric(difftime(DateTime, lag(DateTime), units = "secs")))
  filtered_data <- filtered_data %>% filter(!is.na(Time_Diff))
  if (nrow(filtered_data) < 2) {
    message("Not enough data to calculate burstiness.")
    return(NULL)
  }
  mean_interevent <- mean(filtered_data$Time_Diff)
  sd_interevent <- sd(filtered_data$Time_Diff)
  min_interevent <- min(filtered_data$Time_Diff)
  n <- nrow(filtered_data)
  r <- sd_interevent / mean_interevent
  An_r <- function(r, n) {
    sqrt_n_plus_1 <- sqrt(n + 1)
    sqrt_n_minus_1 <- sqrt(n - 1)
    (sqrt_n_plus_1 * r - sqrt_n_minus_1) / ((sqrt_n_plus_1 - 2) * r + sqrt_n_minus_1)
  }
  novel_burstiness <- An_r(r, n)
  An_y_r <- function(r, n, y) {
    sqrt_n_plus_1 <- sqrt(n + 1)
    sqrt_n_minus_1 <- sqrt(n - 1)
    (n - 2) * (sqrt_n_plus_1 * r - sqrt_n_minus_1 * (1 - n * y)) /
      ((n * sqrt_n_plus_1 - 2 * (n - 1)) * r + sqrt_n_minus_1 * (n - 2 * sqrt_n_plus_1) * (1 - n * y))
  }
  y <- min_interevent / mean(filtered_data$Time_Diff)
  burstiness_min_interevent <- An_y_r(r, n, y)
  return(list(
    Novel_Burstiness = novel_burstiness,
    Burstiness_Min_Interevent = burstiness_min_interevent,
    Number_of_Events = n,
    Mean_Interevent_Time = mean_interevent,
    SD_Interevent_Time = sd_interevent,
    Coefficient_of_Variation = r,
    Min_Interevent_Time = min_interevent
  ))
}

#' Measure Burstiness for All Participants
#'
#' This function calculates burstiness metrics for all participants.
#'
#' @param results The list containing data frames of all participants' data.
#' @return A data frame containing burstiness metrics for all participants.
#' @export
#' @examples
#' burstiness_results <- measure_burstiness(results)
measure_burstiness <- function(results) {
  output_folder <- results$output_folder
  dataframes <- results$dataframes

  burstiness_folder <- file.path(output_folder, "novel_burstiness_measures")
  if (!dir.exists(burstiness_folder)) {
    dir.create(burstiness_folder, recursive = TRUE)
  }

  burstiness_results <- lapply(names(dataframes), function(person) {
    person_dataframes <- dataframes[[person]]
    person_burstiness <- lapply(names(person_dataframes), function(day) {
      df <- person_dataframes[[day]]
      if (is.null(df)) {
        message("Skipping null dataframe.")
        return(NULL)
      }
      burstiness_metrics <- calculate_burstiness(df)
      if (!is.null(burstiness_metrics)) {
        return(data.frame(
          Person = person,
          Day = day,
          Novel_Burstiness = burstiness_metrics$Novel_Burstiness,
          Burstiness_Min_Interevent = burstiness_metrics$Burstiness_Min_Interevent
        ))
      }
    }) %>% bind_rows()

    if (!is.null(person_burstiness) && nrow(person_burstiness) > 0) {
      write.csv(person_burstiness, file.path(burstiness_folder, paste0(person, "_novel_burstiness_metrics.csv")), row.names = FALSE)
    }

    return(person_burstiness)
  }) %>% bind_rows()

  write.csv(burstiness_results, file.path(output_folder, "all_novel_burstiness_metrics.csv"), row.names = FALSE)

  return(burstiness_results)
}


#' Run Full Analysis on GENEActiv Data
#'
#' This function processes all GENEActiv .bin files in the specified folder, computes burstiness metrics for each person,
#' and saves the results in the output folder.
#'
#' @param binfolder The path to the folder containing the .bin files.
#' @param outfolder The path to the output folder where results will be saved.
#' @return A data frame containing burstiness metrics for all subjects.
#' @export
#' @examples
#' burstiness_metrics <- run_full_analysis("D:/Downloads/p004", "D:/Downloads/out/conametest")
run_full_analysis <- function(binfolder, outfolder) {
  results <- compute_for_all(binfolder, outfolder)
  burstiness_metrics <- measure_burstiness(results)
  return(burstiness_metrics)
}
