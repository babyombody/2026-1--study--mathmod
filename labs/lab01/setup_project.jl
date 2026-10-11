#!/usr/bin/env julia
using Pkg
Pkg.add("DrWatson")
using DrWatson
project_name = "project"
initialize_project(project_name; authors="Дмитрий Бахарев", git=false)
println("Проект создан: ", project_name)
println("Перейдите в директорию: cd ", project_name)
