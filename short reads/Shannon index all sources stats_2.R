# ============================================================
# Compare Sample Types SEPARATELY for Time1 and Time2
# ============================================================

cat("\n=== Pairwise Comparisons Between Sample Types (Time1 Only) ===\n")

# Filter for Time1 only
alpha_time1 <- alpha_diversity %>% filter(Time == "Time1")
groups <- unique(alpha_time1$Group)
comparisons <- combn(groups, 2, simplify = FALSE)

results_time1_list <- list()
for(i in seq_along(comparisons)) {
  result <- hutcheson_test(alpha_time1, "Group", 
                           comparisons[[i]][1], comparisons[[i]][2])
  result$TimePoint <- "Time1"
  results_time1_list[[i]] <- result
}

results_time1_df <- do.call(rbind, results_time1_list)
print(results_time1_df)

cat("\n=== Pairwise Comparisons Between Sample Types (Time2 Only) ===\n")

# Filter for Time2 only
alpha_time2 <- alpha_diversity %>% filter(Time == "Time2")
groups <- unique(alpha_time2$Group)
comparisons <- combn(groups, 2, simplify = FALSE)

results_time2_list <- list()
for(i in seq_along(comparisons)) {
  result <- hutcheson_test(alpha_time2, "Group", 
                           comparisons[[i]][1], comparisons[[i]][2])
  result$TimePoint <- "Time2"
  results_time2_list[[i]] <- result
}

results_time2_df <- do.call(rbind, results_time2_list)
print(results_time2_df)

# Combine both for a complete table
cat("\n=== Combined Results ===\n")
all_results <- rbind(results_time1_df, results_time2_df)
print(all_results)