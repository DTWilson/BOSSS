library(hexSticker)

imgurl <- system.file("./man/figures/pixelmech6Big2.png", package="hexSticker")

path_image <- "./man/figures/pixelmech6Big2.png"

library(showtext)
## Loading Google fonts (http://www.google.com/fonts)
font_add_google("Micro 5", db_cache = FALSE)
## Automatically use showtext to render text for future devices
showtext_auto()


sticker(path_image, package="BOSSS", p_size=42, p_y=1.45, p_family = "Micro 5",
        s_x=1, s_y=.85, s_width=.6,
        h_fill="#e0e1dd", h_color="#778da9", p_color = "#a4243b",
        filename="./man/figures/logo.png")
