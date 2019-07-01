g++ -fPIC -O3 -ffast-math -g3 -c -fmessage-length=0 -std=c++14 -MMD -MP -MF"TaylorF2e.d" -MT"TaylorF2e.o" -o "TaylorF2e.o" "TaylorF2e.cpp"
g++ -fPIC -O3 -ffast-math -g3 -c -fmessage-length=0 -std=c++14 -MMD -MP -MF"Amps.d" -MT"Amps.o" -o "Amps.o" "Amps.cpp"
g++ -fPIC -O3 -ffast-math -g3 -c -fmessage-length=0 -std=c++14 -MMD -MP -MF"main.d" -MT"main.o" -o "main.o" "main.cpp"
g++ -fPIC -O3 -ffast-math -g3 -c -fmessage-length=0 -std=c++14 -MMD -MP -MF"wrap.d" -MT"wrap.o" -o "wrap.o" "wrap.cpp"

g++ -shared -fPIC -o "tf2e.so"  ./Amps.o ./TaylorF2e.o ./wrap.o  -lgsl -lgslcblas -lm

g++ -o "Gen_TaylorF2e"  ./Amps.o ./TaylorF2e.o ./main.o  -lgsl -lgslcblas -lm

#run with something like

./Gen_TaylorF2e 8.7 0.25 0.8 0.23 0.23 0.23 0.23 5 5 5 1e18
