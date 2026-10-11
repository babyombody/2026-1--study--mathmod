#!/usr/bin/env julia
# tangle.jl — Генератор отчетов из Literate-скриптов
# Использование: julia tangle.jl <путь_к_скрипту>
#
# Пример: julia --project=. scripts/tangle.jl scripts/01_exponential_growth.jl

using DrWatson
@quickactivate "project"
using Literate

function main()
    if length(ARGS) == 0
        println("""
        Использование: julia tangle.jl <путь_к_скрипту>

        Примеры:
          julia tangle.jl scripts/01_exponential_growth.jl
          julia tangle.jl scripts/02_exponential_growth.jl
        """)
        return
    end

    script_path = ARGS[1]

    if !isfile(script_path)
        error("Файл не найден: $script_path")
    end

    script_name = splitext(basename(script_path))[1]

    println("Генерация из: $script_path")

    # 1. Чистый скрипт (без Markdown-комментариев)
    scripts_dir = scriptsdir(script_name)
    Literate.script(script_path, scripts_dir; credit=false)
    println("  ✓ Чистый скрипт: $(scripts_dir)/$(script_name).jl")

    # 2. Quarto-документ
    quarto_dir = projectdir("markdown", script_name)
    Literate.markdown(
        script_path,
        quarto_dir;
        flavor=Literate.QuartoFlavor(),
        name=script_name,
        credit=false
    )
    println("  ✓ Quarto: $(quarto_dir)/$(script_name).qmd")

    # 3. Jupyter Notebook
    notebooks_dir = projectdir("notebooks", script_name)
    Literate.notebook(
        script_path,
        notebooks_dir;
        name=script_name,
        execute=false,
        credit=false
    )
    println("  ✓ Notebook: $(notebooks_dir)/$(script_name).ipynb")

    println("\nГотово! Все файлы созданы.")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

