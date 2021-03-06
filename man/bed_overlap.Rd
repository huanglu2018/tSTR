% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bed_overlap.R
\name{bed_overlap}
\alias{bed_overlap}
\title{Applied to check the distance between each STR site and closest exon boundaries, if too close, it might be disturbed by alternative splicing(AS), alternative polyadenylation(APA), or alternative promoter usage(APU) events}
\usage{
bed_overlap(query_bed, ref_bed)
}
\arguments{
\item{query_bed}{bed3 or more cols to check their overlap with references}

\item{ref_bed}{bed files as reference, when extracting from gtf, careful about differences between 0-index and 1-index}
}
\value{
combined df with two more columns, "inside_distance_to_ref_boundary" and "outside_distance_to_ref_boundary"; inside_distance_to_ref_boundary=0 means not inside or just at the boundary, outside_distance_to_ref_boundary=0 means not outside or just at the boundary
}
\description{
Applied to check the distance between each STR site and closest exon boundaries, if too close, it might be disturbed by alternative splicing(AS), alternative polyadenylation(APA), or alternative promoter usage(APU) events
}
