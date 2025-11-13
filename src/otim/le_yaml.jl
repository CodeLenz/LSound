#
# Le o arquivo .yaml e retorna valores associados ao 
# problema de otimização
# 
#
function Le_YAML(arquivo::AbstractString,ver=1.0;verbose=false)
    
    # ###########################################################
    # valores padrão (que podem ser modificados via arquivo) 
    # ###########################################################

    # Número de iterações 
    niter = 100

    # Parâmetros do ISLP
    ϵ1 = 0.1
    ϵ2 = 0.1

    # Fração de volume (restrição de volume)
    vf = 0.5

    # Perímetro limite (restrição de perímetro)
    perimetro = 0.0

    # Fator de mudança cheio/vazio
    fatorcv = 5E-3

    # Raio do filtro 
    raio = 0.0

    # Primeiro lemos o arquivo de dados
    dados = YAML.load_file(arquivo)
 
    # Verifica se temos informação sobre a versão do arquivo de dados
    versao = 0.0
    if haskey(dados,"versao")
 
       # Le a versão do arquivo
       versao = dados["versao"]
 
       # Verifica se a versão é compatível
       versao==ver || throw("Le_YAML::versão do arquivo não é compatível com a versão atual") 
         
    end
     
    # Recupera raio do filtro
    if haskey(dados,"raio")

        # recupera como string
        string_raio = dados["raio"]

        # Se foi informado como string, convertemos
        if isa(string_raio,String)
           raio =  parse(Float64,string_raio)
        else
           raio = string_raio
        end

        # Testa consistência da informação 
        raio<0 && throw("Le_YAML::raio deve ser >=0") 

    else
        throw("Le_YAML::raio é uma informação obrigatória") 
    end

    # Recupera ϵ1
    if haskey(dados,"ϵ1")

        # recupera como string
        string_eps1 = dados["ϵ1"]

        # Se foi informado como string, convertemos
        if isa(string_eps1,String)
            ϵ1 =  parse(Float64,string_eps1)
        else
            ϵ1 = string_eps1
        end
 
        # Testa consistência da informação 
        (ϵ1<=0 || ϵ1>=1) && throw("Le_YAML::ϵ1 deve deve estar em (0,1) ") 
        
    else
        println("Parâmetro ϵ1 não foi informado no .yaml. Utilizando o valor padrão ", ϵ1)
    end

    # Recupera ϵ2
    if haskey(dados,"ϵ2")

        # recupera como string
        string_eps2 = dados["ϵ2"]

        # Se foi informado como string, convertemos
        if isa(string_eps2,String)
            ϵ2 =  parse(Float64,string_eps2)
        else
            ϵ2 = string_eps2
        end
 
        # Testa consistência da informação 
        (ϵ2<=0 || ϵ2>=1) && throw("Le_YAML::ϵ2 deve deve estar em (0,1) ") 
        
    else
        println("Parâmetro ϵ2 não foi informado no .yaml. Utilizando o valor padrão ", ϵ2)
    end

    # Recupera o número de iterações
    if haskey(dados,"niter")

        # recupera como string
        string_niter = dados["niter"]

        # Se foi informado como string, convertemos
        if isa(string_niter,String)
            niter =  parse(Int64,string_niter)
        else
            niter = string_niter
        end
 
        # Testa consistência da informação 
        niter>=1 || throw("Le_YAML::Número de iterações deve ser >=1") 
        
    else
        println("Número de iterações não foi informado no .yaml. Utilizando o valor padrão ", niter)
    end

    # Recupera a fração de volume
    if haskey(dados,"volfrac")

        # recupera como string
        string_vf = dados["volfrac"]

        # Se foi informado como string, convertemos
        if isa(string_vf,String)
            vf =  parse(Float64,string_vf)
        else
            vf = string_vf
        end
 
        # Testa consistência da informação 
        (vf<=0||vf>=1) && throw("Le_YAML::Volume fraction deve estar em (0,1) ") 
        
    else
        println("Fração de volume não foi informado no .yaml. Utilizando o valor padrão ", vf)
    end

    # Valor limite para a restrição de perímetro
    if haskey(dados,"perimetro")

        # recupera como string
        string_perimetro = dados["perimetro"]

        # Se foi informado como string, convertemos
        if isa(string_perimetro,String)
            perimetro =  parse(Float64,string_perimetro)
        else
            perimetro = string_perimetro
        end
 
        # Testa consistência da informação 
        (perimetro<0) && throw("Le_YAML::Valor limite do perímetro deve ser maior do que zero") 
        
    else
        println("Perímetro não foi informado no .yaml. Utilizando o valor padrão ", perimetro)
    end


    # Fator de atualização de air/solid
    if haskey(dados,"fatorcv")

        # recupera como string
        string_fatorcv = dados["fatorcv"]

        # Se foi informado como string, convertemos
        if isa(string_fatorcv,String)
            fatorcv =  parse(Float64,string_fatorcv)
        else
            fatorcv = string_fatorcv
        end
 
        # Testa consistência da informação 
        (fatorcv<=0||fatorcv>=1) && throw("Le_YAML::Fator de mudança para cheio/vazio  maior do que zero e menor do que um") 
        
    else
        println("fatorcv não foi informado no .yaml. Utilizando o valor padrão ", fatorcv)
    end

   # Retorna os dados 
   return raio, niter, ϵ1, ϵ2,  vf, perimetro, fatorcv

end