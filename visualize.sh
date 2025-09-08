g++ main.cpp DataParser.cpp -o visualizador -lglfw -lGLEW -lGL -lX11 -lpthread -lXrandr -ldl
./visualizador "./results/60/image_6.00_0.20_0.20_50.00.dat" 6 60
ffmpeg -y -framerate 15 -i output/frame_%05d.png -c:v libx264 -crf 20 -pix_fmt yuv420p animacao_final.mp4 -y
rm output/*