using DGGRID7_jll

function call_dggrid(meta::Dict)
    meta_string = ""
    for (key, val) in meta
        meta_string *= "$(key) $(val)\n"
    end

    tmp_dir = tempname()
    mkdir(tmp_dir)

    meta_path = tmp_dir * "/meta.txt"
    write(meta_path, meta_string)

    DGGRID7_jll.dggrid() do dggrid_path
        run(`$dggrid_path $(meta_path)`)
    end

    # TODO: Remove tmp dir
end