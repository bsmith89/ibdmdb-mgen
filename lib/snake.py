# Here we have a template for aliasing
alias_recipe = "ln -rs {input} {output}"
hardlink_recipe = "ln {input} {output}"
curl_recipe = "curl '{params.url}' > {output}"
curl_unzip_recipe = "curl '{params.url}' | zcat > {output}"
