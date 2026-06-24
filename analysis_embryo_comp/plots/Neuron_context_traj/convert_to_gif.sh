convert \
  -delay 10 \
  -loop 0 \
  -layers Optimize \
  -resize 75% \
  -fuzz 5% \
  -dither FloydSteinberg \
  -colors 64 \
  *.png \
  Reviewer_Figure_1.gif
