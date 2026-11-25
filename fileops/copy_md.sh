#!/bin/bash

SRC="../v62"
NEP="./nep.txt"

# 检查当前目录是否有 nep.txt
if [ ! -f "$NEP" ]; then
    echo "Error: 当前目录没有 nep.txt"
    exit 1
fi

# 复制 v60 中的 md_mv* 文件夹
for src in "$SRC"/md_mv*/; do
    [ -d "$src" ] || continue

    name=$(basename "$src")
    echo "Copying $name ..."

    mkdir -p "$name"

    rsync -a \
        --exclude='nep.txt' \
        --exclude='make_dataset.xyz' \
        --exclude='movie_all.xyz' \
        --exclude='md_mv*.xyz' \
        --exclude='POSCAR_*' \
        --exclude='disp*/' \
        --exclude='struct*/' \
        "$src" "$name/"
done

# 把当前 v62 的 nep.txt 复制到每个计算目录
for md_dir in md_mv*/; do
    [ -d "$md_dir" ] || continue

    for calc_dir in "$md_dir"*/; do
        [ -d "$calc_dir" ] || continue

        cp -f "$NEP" "$calc_dir/nep.txt"
        echo "Copied nep.txt -> $calc_dir"
    done
done

echo "Done."
