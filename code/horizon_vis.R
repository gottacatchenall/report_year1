library(tidyverse)
library(ggthemr)
library(latex2exp)
library(extrafont)
loadfonts()
string_to_tex = function(str) { return(TeX(str)) }

ggthemr('fresh', spacing=3)

forecasting_hor = read.csv('forecasting_horizon.csv')
forecasting_hor %>% 
    group_by(sigma, epsilon) %>%
    summarize(inner=mean(inner_time_until_fail), outer=mean(outer_time_until_fail)) %>%
        ggplot(aes(log(sigma))) + 
        geom_point(aes(y=inner), shape=1, color='blue') +
        geom_point(aes(y=outer), shape=3, color='red') +
         facet_wrap(. ~ epsilon)


df = read.csv('trajectories.csv')
forecasting_hor = df %>% 
    group_by(sigma, epsilon) %>% 
    filter(time %% 20 == T)%>%
    mutate(time = time/100) %>% 
        ggplot(aes(time, innerDifference, color=factor(sigma), shape=factor(sigma))) + 
        geom_line(size=1.5) +
        geom_point(size=3) + 
        theme(aspect.ratio=1.0, 
            text=element_text(family="", size=18),
            axis.title = element_text(size=18),
            axis.text = element_text(family="LM Roman 10", size=16),
            legend.text = element_text(size=16),
            panel.border = element_rect(fill=NA,color="#222222", size=1)
        ) +
        labs(title='Forecasting Horizon for Double Pendulum', x='time', y='Difference between true and forecasted state', color=TeX('Measurement error, $\\sigma$')) +
        scale_color_manual(
            labels=c(expression(10^-12), expression(10^-10), expression(10^-8),expression(10^-6),expression(10^-4),expression(10^-2)),
            values=c("#353945", "#19b4ca", "#40967b", "#ffd678", "#ffae78", "#fb5a3d"),  
        ) +
        scale_shape_manual(values=c(0,1,3,2,5,6))

forecasting_hor

ggsave("forecasting_horizon.png", plot=forecasting_hor, width=12, height =8, device=png())