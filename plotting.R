library(ggplot2)
result <- readr::read_csv("data/processed/boot_all.csv")

result |>
  dplyr::mutate(difference = boot_all - actual) |>
  dplyr::mutate(difference_perc = (boot_all - actual) / actual * 100) |> 
  ggplot(aes(y = difference_perc)) +
  geom_boxplot()

mean(result$boot_all)
mean(result$actual)
