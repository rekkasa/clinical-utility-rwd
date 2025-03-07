library(ggplot2)
result <- readr::read_csv("data/processed/boot_all.csv")
result |>
  dplyr::mutate(difference = boot_all - actual) |>
  dplyr::mutate(
    difference_perc_all = (boot_all - actual) / actual * 100,
    difference_perc_decision = (boot_decision - actual) / actual * 100,
    difference_perc_leftout = (boot_leftout - actual) / actual * 100
  ) |>
  tidyr::pivot_longer(
    cols = tidyr::starts_with("difference_perc"),
    names_to = "difference_type",
    values_to = "difference_value"
  ) |>
  ggplot(aes(y = difference_value, fill = difference_type)) +
  geom_boxplot()

mean(result$boot_all)
mean(result$actual)
