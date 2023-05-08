add_executable(gen_figures gen_figures.cc graphviz.cc tikz.cc)
target_link_libraries(gen_figures satsuma::satsuma)

