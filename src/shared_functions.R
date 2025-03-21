#! /usr/bin/env Rscript

library(tidyverse)
library(magrittr)
library(knitr)
library(kableExtra)
library(cowplot)
library(rlang)

diagnostic_length_plot <- function(data, approx_length=11000, binwidth=200, percentage_y=15000) {
  # Recommended dimensions: fig.width=8, fig.height=3
  # Example call:
  #    diagnostic_length_plot (
  #      data = data,
  #      approx_length = 11000,
  #      binwidth = 200,
  #      percentage_y = 15000
  #     )
  
  perc_length <- c(1.0, 0.90, 0.8, 0.75, 0.70)*approx_length
  max_length <- max(data$length)+1
  
  p <- data %>%
    ggplot2::ggplot(ggplot2::aes(x = length)) +
    ggplot2::geom_histogram(binwidth = binwidth, fill = "skyblue", color = "black") +
    ggplot2::theme_bw() +
    ggplot2::labs(
      x = paste("Nucleotide lengths (binwidth=", binwidth,")", sep = ""),
      y = "Count",
      title = paste("Diagnostic plot to estimate min-length filter for phylogenetic analysis (all data = ", nrow(data), ")", sep = "")
    ) +
    ggplot2::geom_vline(xintercept = perc_length, linetype = "dashed", color = "red") +
    ggplot2::annotate("text", 
             x = perc_length + binwidth*1.5, 
             y = rep(percentage_y, 5),
             label = paste0(round(perc_length), " (", c(100, 90, 80, 75, 70), "%)"),
             hjust = -0.1,
             angle = -90,
             size = 3,
             color = "blue"
    )
  
  return(p)
}

top_x_author_table <- function(data, top_x=25, authors="authors"){
  top_authors <- data %>%
    mutate(
      authors = tolower(!!sym(authors))
    ) %>%
    group_by(authors) %>%
    summarise(
      n = n(),
      regions = paste(unique(region), collapse = ", "),
      countries = paste(unique(country), collapse = ", ")
    ) %>%
    arrange(desc(n)) %>%
    slice_head(n = top_x)
  
  kable(
    top_authors,
    caption = paste("Top ",top_x," most frequent sequence submitters with their region and countries", sep="")
  ) %>%
    kable_styling(
      latex_options = c("hold_position", "striped")
    )
}

top_x_field_table <- function(data, top_x=25, field="host"){
  top_counts <- data %>%
    mutate(
      lower_field = tolower(!!sym(field))
    ) %>%
    group_by(lower_field) %>%
    summarise(
      n = n(),
    ) %>%
    arrange(desc(n)) %>%
    slice_head(n = top_x)
  
  kable(
    top_counts,
    caption = paste0("Top ",top_x," most frequent ", field)
  ) %>%
    kable_styling(
      latex_options = c("hold_position", "striped")
    )
}

diagnostic_time_fill_plot <- function(data, margin_l=0.5, fill="region", title="Frequency and proportion", rel_heights = c(5,7)){
  sdata <- data %>%
    subset(., !is.na(date)) %>%
    subset(., date !="XXXX-XX-XX") %>%
    dplyr::mutate(
      date_adjusted = lubridate::date(gsub("-XX", "-01", date)),
      col_year = substr(date, 1, 4),
      col_mon = substr(date, 6, 7),
      col_day = substr(date, 9, 10)
    ) %>%
    dplyr::group_by(!!sym(fill), col_year) %>%
    dplyr::summarise(
      n=n()
    ) %>%
    dplyr::ungroup(.)
  
  n=sum(sdata$n)
  
  a <- sdata %>%
    ggplot2::ggplot(., ggplot2::aes(x=col_year, y=n, fill=!!sym(fill))) +
    ggplot2::geom_bar(stat="identity", position="stack") +
    ggplot2::theme_bw() + 
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   axis.title.x=ggplot2::element_blank(), axis.title.y=ggplot2::element_blank(),
                   legend.position = "none",
                   plot.margin = margin(l=margin_l,r=0.23,unit="cm"))+
    ggplot2::labs(title=paste(title," (n=", n, ")", sep=""))
  
  b <- sdata %>%
    ggplot2::ggplot(., ggplot2::aes(x=col_year, y=n, fill=!!sym(fill))) +
    ggplot2::geom_bar(stat="identity", position="fill") +
    ggplot2::scale_y_continuous(labels = scales::percent) +
    ggplot2::theme_bw() + 
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle=90, vjust=1, hjust=1, size=6), 
      axis.title.x = ggplot2::element_blank(), axis.title.y=ggplot2::element_blank(),
      legend.position = "bottom", legend.title = ggplot2::element_blank())
  
  cowplot::plot_grid(a, b, ncol=1, rel_heights = rel_heights)
}

diagnostic_time_fill_plot_month <- function(data, margin_l=0.5, fill="region", title="Frequency and proportion", rel_heights = c(5,7)){
  sdata <- data %>%
    subset(., !is.na(date)) %>%
    subset(., date !="XXXX-XX-XX") %>%
    dplyr::mutate(
      date_adjusted = lubridate::date(gsub("-XX", "-01", date)),
      col_year = substr(date, 1, 4),
      col_mon = substr(date, 6, 7),
      col_day = substr(date, 9, 10),
      col_quarter = case_when(
        col_mon %in% c("01", "02", "03") ~ "Q1",
        col_mon %in% c("04", "05", "06") ~ "Q2",
        col_mon %in% c("07", "08", "09") ~ "Q3",
        col_mon %in% c("10", "11", "12") ~ "Q4",
        TRUE ~ col_mon
      ),
      year_mon = paste(col_year, col_mon, sep=":")
    ) %>%
    dplyr::group_by(!!sym(fill), year_mon) %>%
    dplyr::summarise(
      n=n()
    ) %>%
    dplyr::ungroup(.)
  
  n=sum(sdata$n)
  
  a <- sdata %>%
    ggplot2::ggplot(., ggplot2::aes(x=year_mon, y=n, fill=!!sym(fill))) +
    ggplot2::geom_bar(stat="identity", position="stack") +
    ggplot2::theme_bw() + 
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   axis.title.x=ggplot2::element_blank(), axis.title.y=ggplot2::element_blank(),
                   legend.position = "none",
                   plot.margin = margin(l=margin_l,r=0.23,unit="cm"))+
    ggplot2::labs(title=paste(title," (n=", n, ")", sep=""))
  
  b <- sdata %>%
    ggplot2::ggplot(., ggplot2::aes(x=year_mon, y=n, fill=!!sym(fill))) +
    ggplot2::geom_bar(stat="identity", position="fill") +
    ggplot2::scale_y_continuous(labels = scales::percent) +
    ggplot2::theme_bw() + 
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle=90, vjust=1, hjust=1, size=6), 
      axis.title.x = ggplot2::element_blank(), axis.title.y=ggplot2::element_blank(),
      legend.position = "bottom", legend.title = ggplot2::element_blank())
  
  cowplot::plot_grid(a, b, ncol=1, rel_heights = rel_heights)
}


diagnostic_time_fill_plot_quarter <- function(data, margin_l=0.5, fill="region", title="Frequency and proportion", rel_heights = c(5,7)){
  sdata <- data %>%
    subset(., !is.na(date)) %>%
    subset(., date !="XXXX-XX-XX") %>%
    dplyr::mutate(
      date_adjusted = lubridate::date(gsub("-XX", "-01", date)),
      col_year = substr(date, 1, 4),
      col_mon = substr(date, 6, 7),
      col_day = substr(date, 9, 10),
      col_quarter = case_when(
        col_mon %in% c("01", "02", "03") ~ "Q1",
        col_mon %in% c("04", "05", "06") ~ "Q2",
        col_mon %in% c("07", "08", "09") ~ "Q3",
        col_mon %in% c("10", "11", "12") ~ "Q4",
        TRUE ~ col_mon
      ),
      year_mon = paste(col_year, col_quarter, sep=":")
    ) %>%
    dplyr::group_by(!!sym(fill), year_mon) %>%
    dplyr::summarise(
      n=n()
    ) %>%
    dplyr::ungroup(.)
  
  n=sum(sdata$n)
  
  a <- sdata %>%
    ggplot2::ggplot(., ggplot2::aes(x=year_mon, y=n, fill=!!sym(fill))) +
    ggplot2::geom_bar(stat="identity", position="stack") +
    ggplot2::theme_bw() + 
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   axis.title.x=ggplot2::element_blank(), axis.title.y=ggplot2::element_blank(),
                   legend.position = "none",
                   plot.margin = margin(l=margin_l,r=0.23,unit="cm"))+
    ggplot2::labs(title=paste(title," (n=", n, ")", sep=""))
  
  b <- sdata %>%
    ggplot2::ggplot(., ggplot2::aes(x=year_mon, y=n, fill=!!sym(fill))) +
    ggplot2::geom_bar(stat="identity", position="fill") +
    ggplot2::scale_y_continuous(labels = scales::percent) +
    ggplot2::theme_bw() + 
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle=90, vjust=1, hjust=1, size=6), 
      axis.title.x = ggplot2::element_blank(), axis.title.y=ggplot2::element_blank(),
      legend.position = "bottom", legend.title = ggplot2::element_blank())
  
  cowplot::plot_grid(a, b, ncol=1, rel_heights = rel_heights)
}

earliest_records_table <- function(data, first=5, date="date", other_fields="accession region"){
  other_cols=strsplit(other_fields, " ")[[1]]
  
  earliest_table <- data %>%
    filter(!is.na(!!sym(date))) %>%
    filter(!!sym(date) != "XXXX-XX-XX") %>%
    mutate(
      date_adjusted = lubridate::date(gsub("-XX", "-01", !!sym(date))),
      col_year = substr(!!sym(date), 1, 4),
      col_mon = substr(!!sym(date), 6, 7),
      col_day = substr(!!sym(date), 9, 10)
    ) %>%
    arrange(date_adjusted) %>%
    select(date_adjusted, !!!syms(other_cols)) %>%
    slice_head(n = first)
  
  kable(
    earliest_table,
    caption = paste0("Top ", first, " earliest records")
  ) %>%
    kable_styling(
      latex_options = c("hold_position", "striped")
    )
}

