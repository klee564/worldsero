
train_list <- train_input(1,del=21)
inputlist <- input_multinom()
mydat <- preprocVE2()  %>% df_to_standat2


library(future)
library(furrr)
library(magrittr)
plan(multisession, workers = 10)


simul_list <- furrr::future_map(1:100,~runsimul(train_list,inputlist,mydat),
                               .options = furrr_options(seed = 123))

simul_list %>% do.call("rbind",.) %>% group_by(param) %>%
  summarise(mean(coverage)) %>% View()
simul_list %>% do.call("rbind",.) %>% dplyr::filter(param %in% paste0("mu",1:7)) %>%
  pull(coverage) %>% mean
simul_list %>% do.call("rbind",.) %>% dplyr::filter(param %in% paste0("pmu",1:7)) %>%
  pull(coverage) %>% mean
simul_list %>% do.call("rbind",.) %>% dplyr::filter(param %in% c(paste0("mu",1:7),paste0("pmu",1:7))) %>%
  pull(coverage) %>% mean
simul_list %>% do.call("rbind",.) %>% dplyr::filter(param %in% c("alpha","beta","beta_multi")) %>%
  group_by(param) %>% summarise(mean(coverage))

saveRDS(simul_list,"simul_list3.rds")
simul_list %>% do.call("rbind",.) %>% group_by(param) %>%
  summarise(mean(err/truval^2))

simul_list %>% do.call("rbind",.) %>% dplyr::filter(param=="mu") %>% head

as.data.frame(table(mydat$g2))


