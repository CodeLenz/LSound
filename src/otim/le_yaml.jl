#
# Le o arquivo .yaml e retorna valores associados ao 
# problema de otimização
# 
#
function Le_YAML(arquivo::AbstractString,ver=1.0;verbose=false)
    
    #
    # valores padrão (que podem ser modificados via arquivo) 
    #

    # Evolution Rate
    er = 0.01

    # Número de iterações 
    niter = 100

    # Número de iterações para o histórico de sensibilidade
    nhisto = 1

    # Valor padrão de parametrização 
    # PEREIRA ou DUHRING
    param = "PEREIRA"

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
    raio = 0.0
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

    # Recupera Evolution Rate (er)
    if haskey(dados,"er")

        # recupera como string
        string_er = dados["er"]

        # Se foi informado como string, convertemos
        if isa(string_er,String)
            er =  parse(Float64,string_er)
        else
            er = string_er
        end
 
        # Testa consistência da informação 
        (er<=0||er>=1) && throw("Le_YAML::Evolution Rate raio deve estar em (0,1) ") 
        
    else
        println("Evolution Rate não foi informado no .yaml. Utilizando o valor padrão ", er)
    end

     # Recupera o número de iterações
     if haskey(dados,"niter")

        # recupera como string
        string_niter = dados["niter"]

        # Se foi informado como string, convertemos
        if isa(string_er,String)
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
        println("Evolution Rate não foi informado no .yaml. Utilizando o valor padrão ", er)
    end

    # Recupera a fração de volume
    if haskey(dados,"parametrizacao")

        # recupera como string
        parametrizacao = dados["parametrizacao"]

        # verifica se é uma das parametrizações válidas
        parametrizacao in ["PEREIRA"; "DUHRING"]  || error("Parametrização $parametrizacao é inválida") 

    else
        println("parametrização não foi informada no .yaml. Utilizando o valor padrão ", parametrizacao)
    end


    # Recupera o número de iterações para o histórico de sensibilidades
     if haskey(dados,"nhisto")

        # recupera como string
        string_nhisto = dados["nhisto"]

        # Se foi informado como string, convertemos
        if isa(string_nhisto,String)
            nhisto =  parse(Int64,string_nhisto)
        else
            nhisto = string_nhisto
        end
 
        # Testa consistência da informação 
        nhisto>=1 || throw("Le_YAML::Número de iterações para o histórico de sensibilidades deve ser >=1") 
        
    else
        println("Número de iterações para o histórico de sensibilidade não foi informado no .yaml. Utilizando o valor padrão ", nhisto)
    end


   # Retorna os dados 
   return raio, niter, nhisto, er, vf, parametrizacao

end