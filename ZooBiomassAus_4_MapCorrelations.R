

Net <- readRDS("Biomass_Layer_Nets.rds")
Net <- calc(Net, fun = mean)
names(Net) <- "Nets"
Net_df <- as.data.frame(rasterToPoints(Net))

CPR <- readRDS("Biomass_Layer_CPR.rds")
CPR <- calc(CPR, fun = mean)
names(CPR) <- "CPR"
CPR_df <- as.data.frame(rasterToPoints(CPR))

LOPC <- readRDS("Biomass_Layer_LOPC.rds")
LOPC <- calc(LOPC, fun = mean)
names(LOPC) <- "LOPC"
LOPC_df <- as.data.frame(rasterToPoints(LOPC))

All <- readRDS("Biomass_Layer_Joint.rds")

df <- Net_df

df <- left_join(df, LOPC_df, by = c("x", "y"))
df <- left_join(df, CPR_df, by = c("x", "y"))


ggplot(df, aes(Nets, CPR)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  xlab("Nets") +
  ylab("CPR") +
  theme_bw() +
  theme(aspect.ratio=1,
        axis.title = element_text(size=12),
        axis.text = element_text(size=12))
mdl <- lm(Nets ~ CPR, data = df)
summary(mdl)
ggsave("Figures/Corr_NetsCPR.png", dpi=300, scale=1)


ggplot(df, aes(Nets, LOPC)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  xlab("Nets") +
  ylab("LOPC") +
  theme_bw() +
  theme(aspect.ratio=1,
        axis.title = element_text(size=12),
        axis.text = element_text(size=12))
ggsave("Figures/Corr_NetsLOPC.png", dpi=300, scale=1)

mdl <- lm(Nets ~ LOPC, data = df)



ggplot(df, aes(CPR, LOPC)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  xlab("CPR") +
  ylab("LOPC") +
  theme_bw() +
  theme(aspect.ratio=1,
        axis.title = element_text(size=12),
        axis.text = element_text(size=12))
ggsave("Figures/Corr_CPRLOPC.png", dpi=300, scale=1)

mdl <- lm(CPR ~ LOPC, data = df)

