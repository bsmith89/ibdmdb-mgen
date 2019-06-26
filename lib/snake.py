# Here we have a template for aliasing
alias_recipe = "ln -rs {input} {output}"
alias_fmt = lambda input, output: alias_recipe.format(input=input, output=output)
curl_recipe = "sleep 1; curl '{params.url}' > {output}"
curl_unzip_recipe = "sleep 1; curl '{params.url}' | zcat > {output}"
