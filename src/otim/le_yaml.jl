#
# Le o arquivo .yaml e retorna valores associados ao 
# problema de otimização
# 
#
function Le_YAML(arquivo::AbstractString,ver=1.0;verbose=false)
    
    #
    # valores padrão (que podem ser modificados via arquivo) 
    #

    # Número de iterações 
    niter = 100

    # Número de iterações para o histórico de sensibilidade
    nhisto = 1

    # Parâmetros do ISLP
    ϵ1 = 0.1
    ϵ2 = 0.1

    # Valor mínimo da variável de projeto 
    γmin = 0.0

    # Valor máximo da variável de projeto 
    γmax = 1

    # Valor padrão de parametrização 
    # PEREIRA ou DUHRING
    parametrizacao = "DUHRING"

    # Valor do 'ponto de partida', lembrando que não podemos começar com todas as posições nulas, pois 
    # isso vai fazer com que a atualização de volume seja  0*(1+er) = sempre zero. Então, podemos começar 
    # com um padrão que seja fisicamente adequado para o problema em questão.
    partida = 1.0

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
        if isa(string_er,String)
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

    # Recupera o valor mínimo da variável de projeto
    if haskey(dados,"γmin")

        # recupera como string
        string_γmin = dados["γmin"]

        # Se foi informado como string, convertemos
        if isa(string_γmin,String)
            γmin =  parse(Float64,string_γmin)
        else
            γmin = string_γmin
        end
 
        # Testa consistência da informação 
        ( 0<= γmin <=1 ) || throw("Le_YAML::γmin deve estar em (0,1) ") 
        
    else
        println("γmin não foi informado no .yaml. Utilizando o valor padrão ", γmin)
    end


    # Recupera o valor máximo da variável de projeto
    if haskey(dados,"γmax")

        # recupera como string
        string_γmax = dados["γmax"]

        # Se foi informado como string, convertemos
        if isa(string_γmax,String)
            γmax =  parse(Float64,string_γmax)
        else
            γmax = string_γmax
        end
 
        # Testa consistência da informação 
        (0<γmax<=1) || throw("Le_YAML::γmax deve estar em (0,1) ") 
        
    else
        println("γmax não foi informado no .yaml. Utilizando o valor padrão ", γmax)
    end

    # Teste para ver se γmin < γmax
    (γmin<γmax) || throw("Le_YAML::γmin deve ser menor do que γmax") 

    # Recupera o tipo de parametrização do material
    if haskey(dados,"parametrizacao")

        # recupera como string
        parametrizacao = dados["parametrizacao"]

        # verifica se é uma das parametrizações válidas
        parametrizacao in ["PEREIRA"; "DUHRING"]  || error("Parametrização $parametrizacao é inválida") 

    else
        println("Parametrização não foi informada no .yaml. Utilizando o valor padrão ", parametrizacao)
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

    # Recupera o valor do 'ponto de partida'
    if haskey(dados,"partida")

        # recupera como string
        string_partida = dados["partida"]

        # Se foi informado como string, convertemos
        if isa(string_partida,String)
            partida =  parse(Float64,string_partida)
        else
            partida = string_partida
        end
 
        # Testa consistência da informação 
        ( 0<= partida <=1) || throw("Le_YAML::partida deve estar em (0,1) ") 
        
    else
        println("Ponto de partida não foi informado no .yaml. Utilizando o valor padrão ", partida)
    end



   # Retorna os dados 
   return raio, niter, nhisto, ϵ1, ϵ2,  vf, parametrizacao, γmin, γmax, partida 

end