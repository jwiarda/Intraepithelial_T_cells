library(vegan)
library(funfuns)
library(tidyverse)

setwd('~/Jayne/')

set.seed(42)
CD3 <- read_csv('new_plots/NG8_FC_control_NMDS_CD27all.csv')

CD3$Time <- ifelse(CD3$Time == 7, 4, 
                   ifelse(CD3$Time == 21, 6, 
                          ifelse(CD3$Time == 35, 8, NA)))
CD3$Time <- factor(CD3$Time, levels = c('4','6','8'))
CD3$pignum <- substr(CD3$Sample, 0, 3)
CD3$Sample_ID <- paste(CD3$pignum, CD3$Tissue, CD3$Time, sep = "_")
CD3$set <- paste(CD3$Tissue, CD3$Time, sep = "_")
CD3$Tissue <- factor(CD3$Tissue, levels = c('Jejunum','Ileum','Cecum','Colon'))# set the order of the tissues

##

CD3_gath <- CD3 %>% select(-Sample, -Diet, -set, -Sample_ID, -pignum) %>%
  gather(key='cell_type', value='counts', -c(Tissue, Time))

CD3_sum <- CD3_gath %>%
  mutate(cell_type = sub('Q[1-4]: ', '', cell_type), 
         cell_type = sub(' , ', '', cell_type),
         cell_type = sub('pos', '\\+', cell_type),
         cell_type = gsub('/', '', cell_type),
         counts=as.numeric(counts)) %>%
  group_by(cell_type, Tissue) %>%
  summarise(total_counts = sum(counts)) %>% 
  group_by(Tissue) %>%
  mutate(percent = total_counts/sum(total_counts) * 100) %>% 
  select(-total_counts) %>% spread(key=Tissue, value = percent)

CD3_sum %>% write_csv('Jayne_flow_cell_types.csv')

##

rownames(CD3) <- CD3$Sample_ID  # sets rownames

CD3_meta <- CD3 %>% select(-(5:20)) # meta data
CD3_data <- CD3 %>% select(5:20) # population frequency data

# now we will express everything as a percent of CD3+ events
CD3_data <- (CD3_data / rowSums(CD3_data)) * 100

rowSums(CD3_data) # should all be equal to 100

NMDS <- NMDS_ellipse(metadata = CD3_meta, OTU_table = CD3_data, grouping_set = 'set')

stressplot(NMDS[[3]])
NMDS[[3]]$stress # check stress level according to parameters below:
# stress > 0.2 = bad
# stress < 0.1 =  good (ordination distances do a good job at representing underlying distances)

# extract some metadata from the set variable in the ellipse dataframe
NMDS[[2]] <- NMDS[[2]] %>%
  mutate(Tissue = sub('(.*)_(.*)','\\1',group), 
         Time = sub('(.*)_(.*)','\\2',group)) %>%
  mutate(Tissue = factor(Tissue, levels = c('Jejunum','Ileum','Cecum','Colon')), 
         Time = factor(Time, levels = c('4','6','8'))) # set the order of the tissues


#F8766D cecum
#5BB300 colon
#C77CFF Ileum
#00BFC4 Jejunum

cols <- c('#00BFC4','#C77CFF','#F8766D','#5BB300')

# cols <- c("#F8766D", "#F8766D","#F8766D","#5BB300", "#5BB300","#5BB300", "#C77CFF" ,"#C77CFF","#C77CFF", "#00BFC4", "#00BFC4", "#00BFC4")

NMDS[[1]] %>% arrange(Time) %>% 
  ggplot(aes(x=MDS1, y=MDS2, fill=Tissue, colour = Tissue)) +
  geom_polygon(data = NMDS[[2]], aes(x=NMDS1, y=NMDS2, fill=Tissue, group=group), inherit.aes = FALSE, alpha=0.2)+
  geom_point(shape =21, size=0.75, alpha = 0.2) + 
  geom_point(aes(x=centroidX, y=centroidY, shape=Time), size=3)+
  geom_path(aes(x=centroidX, y=centroidY, group=Tissue))+
  geom_segment(aes(xend=centroidX, yend=centroidY, color=Tissue), linejoin = 'round', size=0.5, alpha = 0.2)+
  ggtitle('NMDS ordination of T-cell community similarities', subtitle = 'stress = 0.099')+
  scale_fill_manual(values=c(cols), aesthetics = c('fill', 'color'))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))


##########

# PERMANOVA #
# do certain factors influence CD3+ community similarity?
# CD3_meta$Time <- as.character(CD3_meta$Time)
adon <- adonis(formula = CD3_data ~ Tissue * Time, data = CD3_meta)
adon

#

pwadon <- pairwise.adonis(CD3_data, CD3_meta$set, p.adjust.m = 'none')

pwadon_time <- pwadon[grep('(.*)_(.*) vs (\\1)_(.*)', pwadon$pairs),]
pwadon_tissue <- pwadon[grep('(.*)_(.*) vs (.*)_(\\2)', pwadon$pairs),]

all_pw <- rbind(pwadon_time, pwadon_tissue)

all_pw$p.adjusted <- p.adjust(all_pw$p.value, method = 'fdr')
all_pw$p.plot <- signif(all_pw$p.adjusted, 3)

all_pw

#################################################


####### insertion ######

set.seed(42)
CD3 <- read_csv('new_plots/NG8_FC_control_NMDS_8subsets.csv')

CD3$Time <- ifelse(CD3$Time == 7, 4, 
                   ifelse(CD3$Time == 21, 6, 
                          ifelse(CD3$Time == 35, 8, NA)))
CD3$Time <- factor(CD3$Time, levels = c('4','6','8'))
CD3$pignum <- substr(CD3$Sample, 0, 3)
CD3$Sample_ID <- paste(CD3$pignum, CD3$Tissue, CD3$Time, sep = "_")
CD3$set <- paste(CD3$Tissue, CD3$Time, sep = "_")
CD3$Tissue <- factor(CD3$Tissue, levels = c('Jejunum','Ileum','Cecum','Colon'))# set the order of the tissues
# CD3$Time <- factor(CD3$Time, levels = c('7','21','35'))



rownames(CD3) <- CD3$Sample_ID  # sets rownames

CD3_meta <- CD3 %>% select(-(5:12)) # meta data
CD3_data <- CD3 %>% select(5:12) # population frequency data

# now we will express everything as a percent of CD3+ events
CD3_data <- (CD3_data / rowSums(CD3_data)) * 100

rowSums(CD3_data) # should all be equal to 100

NMDS <- NMDS_ellipse(metadata = CD3_meta, OTU_table = CD3_data, grouping_set = 'set')

# NMDS_ellipse() is a wrapper function i wrote that makes plotting these ordinations easy in ggplot
# it uses 

stressplot(NMDS[[3]])
NMDS[[3]]$stress # check stress level according to parameters below:
# stress > 0.2 = bad
# stress < 0.1 = really good (ordination distances do a good job at representing underlying distances)

# extract some metadata from the set variable in the ellipse dataframe
NMDS[[2]] <- NMDS[[2]] %>%
  mutate(Tissue = sub('(.*)_(.*)','\\1',group), 
         Time = sub('(.*)_(.*)','\\2',group)) %>%
  mutate(Tissue = factor(Tissue, levels = c('Jejunum','Ileum','Cecum','Colon')), 
         Time = factor(Time, levels = c('4','6','8'))) # set the order of the tissues


#F8766D cecum
#5BB300 colon
#C77CFF Ileum
#00BFC4 Jejunum

cols <- c('#00BFC4','#C77CFF','#F8766D','#5BB300')

# cols <- c("#F8766D", "#F8766D","#F8766D","#5BB300", "#5BB300","#5BB300", "#C77CFF" ,"#C77CFF","#C77CFF", "#00BFC4", "#00BFC4", "#00BFC4")

NMDS[[1]] %>% arrange(Time) %>% 
  ggplot(aes(x=MDS1, y=MDS2, fill=Tissue, colour = Tissue)) +
  geom_polygon(data = NMDS[[2]], aes(x=NMDS1, y=NMDS2, fill=Tissue, group=group), inherit.aes = FALSE, alpha=0.2)+
  geom_point(shape =21, size=0.75, alpha = 0.2) + 
  geom_point(aes(x=centroidX, y=centroidY, shape=Time), size=3)+
  geom_path(aes(x=centroidX, y=centroidY, group=Tissue))+
  geom_segment(aes(xend=centroidX, yend=centroidY, color=Tissue), linejoin = 'round', size=0.5, alpha = 0.2)+
  ggtitle('NMDS ordination of T-cell community similarities', subtitle = 'stress = 0.065')+
  scale_fill_manual(values=c(cols), aesthetics = c('fill', 'color'))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

##########

# PERMANOVA #
# do certain factors influence CD3+ community similarity?
CD3_meta$Time <- as.character(CD3_meta$Time)
adon <- adonis(formula = CD3_data ~ Tissue * Time, data = CD3_meta)
adon

#

pwadon <- pairwise.adonis(CD3_data, CD3_meta$set, p.adjust.m = 'none')

pwadon_time <- pwadon[grep('(.*)_(.*) vs (\\1)_(.*)', pwadon$pairs),]
pwadon_tissue <- pwadon[grep('(.*)_(.*) vs (.*)_(\\2)', pwadon$pairs),]

all_pw <- rbind(pwadon_time, pwadon_tissue)

all_pw$p.adjusted <- p.adjust(all_pw$p.value, method = 'fdr')
all_pw$p.plot <- signif(all_pw$p.adjusted, 3)

all_pw
