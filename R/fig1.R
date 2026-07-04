#bar_data <- data.frame(
#  Category = factor(
#    c("MS", "MS/MS", "SIRIUS",
#      "MS", "MS/MS", "SIRIUS"),
#    levels = c("MS", "MS/MS", "SIRIUS")
#  ),
#  Mode = factor(
#    c("Negative", "Negative", "Negative",
#      "Positive", "Positive", "Positive"),
#    levels = c("Negative", "Positive")
#  ),
#  Features = c(14756, 10268, 7855,
#               25587, 10484, 7919)
#)
#
#feature_summary_plot(
#    bar_data,
#    outfile = "Figure_features_summary2.tiff"
#
#)

feature_summary_plot <- function(
    bar_data,
    outfile = NULL,
    width = 7.2,
    height = 6.8,
    dpi = 600,
    neg_col = "#C65D3A",
    pos_col = "#465C88",
    neu_col = "#D9D9D9",
    txt_col = "#1A1A1A",
    subtle_col = "#6E6E6E",
    bg_col = "white"
){

    ##------------------------------------------------------------
    ## 自动计算Y轴
    ##------------------------------------------------------------

    ymax <- max(bar_data$Features, na.rm = TRUE)

    upper <- pretty(c(0, ymax * 1.08), n = 3)
    ylim_max <- max(upper)
    upper <- upper[upper <= ymax]

    ##------------------------------------------------------------
    ## 数据整理
    ##------------------------------------------------------------

    bar_data <- bar_data |>
        dplyr::group_by(Mode) |>
        dplyr::mutate(
            Percent_vs_MS =
                Features / Features[Category == "MS"] * 100,
            Label = scales::comma(Features)
        ) |>
        dplyr::ungroup()

    make_donut_data <- function(bar_data, mode){

        tmp <- bar_data[bar_data$Mode == mode, ]

        msms <- tmp$Features[tmp$Category == "MS/MS"]
        sirius <- tmp$Features[tmp$Category == "SIRIUS"]

        rate <- if (length(msms) == 0 || msms == 0) {
            0
        } else {
            round(sirius / msms * 100, 1)
        }

        data.frame(
            Group = c("Annotated", "Unannotated"),
            Value = c(rate, 100 - rate)
        )
    }


    p_bar <-
        ggplot2::ggplot(
            bar_data,
            ggplot2::aes(
                x = Category,
                y = Features,
                fill = Mode
            )
        ) +
        ggplot2::geom_col(
            position = ggplot2::position_dodge(width = 0.7),
            width = 0.65,
            color = NA
        ) +
        ggplot2::geom_text(
            ggplot2::aes(label = scales::comma(Features)),
            position = ggplot2::position_dodge(width = 0.72),
            vjust = -0.35,
            size = 4.5,
            fontface = "bold",
            colour = txt_col
        ) +
        ggplot2::scale_fill_manual(
            values = c(
                Negative = neg_col,
                Positive = pos_col
            )
        ) +
        ggplot2::scale_y_continuous(
            breaks = upper,
            labels = scales::comma,
            expand = ggplot2::expansion(mult = c(0, 0.04))
        ) +
        ggplot2::coord_cartesian(
            ylim = c(0, ylim_max),
            clip = "off"
        ) +
        ggplot2::labs(
            x = NULL,
            y = "Number of features",
            fill = NULL
        ) +
        ggplot2::theme_minimal(base_size = 13) +
        ggplot2::theme(

            legend.position = c(0.17,0.97),
            legend.justification = c(0,1),
            legend.direction = "horizontal",
            legend.text = ggplot2::element_text(
                size = 12,
                colour = txt_col
            ),
            legend.key.width = grid::unit(1.1,"cm"),
            legend.key.height = grid::unit(0.45,"cm"),
            legend.background = ggplot2::element_blank(),

          axis.title.y = ggplot2::element_text(
              size = 14,
              face = "bold",
              colour = txt_col
          ),

          axis.text.x = ggplot2::element_text(
              size = 12.5,
              face = "bold",
              colour = txt_col
          ),

          axis.text.y = ggplot2::element_text(
              size = 11.5,
              colour = subtle_col
          ),

          panel.grid.major.x = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(),

          panel.grid.major.y =
              ggplot2::element_line(
                  colour="#ECECEC",
                  linewidth=0.5
              ),

          axis.line.x =
              ggplot2::element_line(
                  colour=txt_col,
                  linewidth=0.5
              ),

          axis.line.y = ggplot2::element_blank(),

          axis.ticks.x =
              ggplot2::element_line(
                  colour=txt_col,
                  linewidth=0.4
              ),

          axis.ticks.y = ggplot2::element_blank(),

          plot.background =
              ggplot2::element_rect(
                  fill=bg_col,
                  colour=NA
              ),

          panel.background =
              ggplot2::element_rect(
                  fill=bg_col,
                  colour=NA
              ),

          plot.margin =
              ggplot2::margin(
                  25,18,12,12
              )
      )

    ##------------------------------------------------------------
    ## donut函数
    ##------------------------------------------------------------

    make_donut_only <- function(df, main_col){

        df <- df |>
            dplyr::mutate(

                Group = factor(
                    Group,
                    levels = c(
                        "Annotated",
                        "Unannotated"
                    )
                ),

                fraction = Value/sum(Value),

                ymax = cumsum(fraction),

                ymin = c(0,head(ymax,-1))
            )

        ggplot2::ggplot(
            df,
            ggplot2::aes(
                ymax=ymax,
                ymin=ymin,
                xmax=4,
                xmin=3,
                fill=Group
            )
        )+
            ggplot2::geom_rect(
                colour="white",
                linewidth=1.2
            )+
            ggplot2::scale_fill_manual(
                values=c(
                    Annotated=main_col,
                    Unannotated=neu_col
                )
            )+
            ggplot2::coord_polar(
                theta="y",
                clip="off"
            )+
            ggplot2::xlim(c(1.35,4.25))+
            ggplot2::theme_void()+
            ggplot2::theme(
                legend.position="none",
                plot.background=
                    ggplot2::element_rect(fill=NA,colour=NA),
                panel.background=
                    ggplot2::element_rect(fill=NA,colour=NA)
            )

    }

    make_donut <- function(df, colour){

        pct <- round(df$Value[df$Group=="Annotated"],1)

        ring <- make_donut_only(df, colour)

        cowplot::ggdraw(ring)+
            cowplot::draw_label(
                paste0(pct,"%"),
                x=.5,
                y=.5,
                fontface="bold",
                size=14,
                colour=txt_col
            )+
            cowplot::draw_label(
                "Annotation rate",
                x=.5,
                y=0,
                fontface="bold",
                size=11,
                colour=subtle_col
            )
    }

    donut_neg <- make_donut_data(bar_data, "Negative")
    donut_pos <- make_donut_data(bar_data, "Positive")

    p_neg <- make_donut(donut_neg, neg_col)
    p_pos <- make_donut(donut_pos, pos_col)

    final_plot <-
        cowplot::ggdraw()+
        cowplot::draw_plot(
            p_bar,
            0,0,1,1
        )+
        cowplot::draw_plot(
            p_neg,
            .46,.65,.24,.24
        )+
        cowplot::draw_plot(
            p_pos,
            .68,.65,.24,.24
        )

    if(!is.null(outfile)){

        ggplot2::ggsave(
            filename = outfile,
            plot = final_plot,
            width = width,
            height = height,
            dpi = dpi,
            compression = "lzw",
            bg = "white"
        )

    }

    return(final_plot)

}


