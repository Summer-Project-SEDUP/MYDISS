% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hello.R
\name{compute_one_person}
\alias{compute_one_person}
\title{Compute Data for One Participant}
\usage{
compute_one_person(file_name, binfolder, outfolder)
}
\arguments{
\item{file_name}{The name of the .bin file for the participant.}

\item{binfolder}{The path to the folder containing the .bin files.}

\item{outfolder}{The path to the output folder where results will be saved.}
}
\value{
A list containing data frames for each day for that participant.
}
\description{
This function processes a single .bin file for a participant, computes posture changes, and saves results.
}
\examples{
person_dataframes <- compute_one_person("participant.bin", "path/to/binfolder", "path/to/outfolder")
}
