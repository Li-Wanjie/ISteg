function [IA_stego,trace] = IA_Post(cover_image,stego_image)
%IA_POST 免疫算法对隐写图像进行免疫处理操作
%   自适应得到隐写后处理之后较好的图像

    tic;
    cover = double(cover_image);
    stego = double(stego_image);
    [M, N] = size(stego);
    imgRes = stego - cover;                     % 载体和载密之间的的残差
    stego = reshape(stego, [M*N,1]);
    imgRes = reshape(imgRes,[M*N,1]);
    modifyRange = find(imgRes == 1 | imgRes == -1);      % 取载密图像修改点的索引值

    Anti_D = M * N;                             % 免疫个体维数,M和N是图像尺寸大小
    Population_Num = 30;                       % 免疫个体数目(种群中抗体数量)
    Anti_Gengs_Upper = 255;                     % 抗体变异后取值的上界
    Anti_Gengs_Lower = 0;                       % 抗体变异后取值的下界
    IA_Maximum = 100;                           % 最大免疫代数,一般取值范围100-500
    gen = 1;                                    % 免疫代数
    Ncl = 40;                                   % 克隆个数
    C = 4;                                      %  隐写后处理单元加减4，不影响秘密信息提取
    randNumtemp = 0;

    %% 初始化种群
    f = repmat(stego, 1, Population_Num);
    cond_anti = cell(1,Population_Num);
    Fitness = zeros(1,Population_Num);

    %% 计算抗体-抗原亲和度
    for Anti = 1:Population_Num    
        cond_anti(Anti) = {reshape(f(:,Anti),M,N)};
        Fitness(Anti) = Distance(cover, cond_anti{Anti});    
    end
    % 抗体-抗原亲和度按升序排列
    [SortFitness, Index]=sort(Fitness);
    Sortf=f(:, Index);

    %% 免疫循环
    while gen <= IA_Maximum
        for i = 1:Population_Num/10
            %%%%%选抗体-抗原亲和度前Population_Num/10个体进行免疫操作%%%%%%
            antibody = Sortf(:,i);
            clone_antibody = repmat(antibody,1,Ncl);        % 抗体克隆，Ncl=10
            for j = 1:Ncl 
                %%%%%%%%%%%%%%%%%变异%%%%%%%%%%%%%%%%%%%
                randNum(gen) = modifyRange(randi(numel(modifyRange),1,1));      % 取随机数，表示随机变异抗体的1个单元
                if randNum(gen) ~= randNumtemp  % 判断是否与上一个变异像素点位置相同
                    alfa = randsrc(1,1,[-1 1]); % 产生-1或1，随机对嵌入的单元进行加减
                    if (clone_antibody(randNum(gen),j) + alfa * C) >= Anti_Gengs_Lower & (clone_antibody(randNum(gen),j) + alfa * C) <= Anti_Gengs_Upper
                        clone_antibody(randNum(gen),j) = clone_antibody(randNum(gen),j) + alfa * C;
                    end
                end
                randNumtemp = randNum(gen);
            end
            clone_antibody(:,1) = Sortf(:,i);               %保留克隆源个体
            %%%%%%%%%%克隆抑制，保留亲和度最高的个体%%%%%%%%%%
            cond_clone_antibody = cell(1,Ncl);
            clone_antibodyFitness = zeros(1,Ncl);
            for j = 1:Ncl     
                cond_clone_antibody(j) = {reshape(clone_antibody(:,j),M,N)};
                clone_antibodyFitness(j) = Distance(cover,cond_clone_antibody{j});
            end
            [SortCloneAntiFitness,Index]=sort(clone_antibodyFitness);
            aMSLL(i)=SortCloneAntiFitness(1);
            clone_antibodySortf=clone_antibody(:,Index);
            af(:,i)=clone_antibodySortf(:,1);
        end 
        %%%%%%%%%%%%%%%%%%%%%%%种群刷新%%%%%%%%%%%%%%%%%%%%%%%
        bf = repmat(stego, 1, 9*(Population_Num/10));   % 新抗体加入
        for w = 1:9*(Population_Num/10) 
            %%%%%%%%%%%%%%%%%随机抗体%%%%%%%%%%%%%%%%%%%
            randNum1 = modifyRange(randperm(numel(modifyRange),gen));      % 取随机数，表示随机变异抗体的gen个单元
            for q = 1:gen
                alfa1 = randsrc(1,1,[-1 1]);           % 产生-1或1，随机对嵌入的单元进行加减
                if (bf(randNum1(q),w) + alfa1 * C) >= Anti_Gengs_Lower & (bf(randNum1(q),w) + alfa1 * C) <= Anti_Gengs_Upper
                    bf(randNum1(q),w) = bf(randNum1(q),w) + alfa1 * C;
                end
            end
        end
        new_anti = cell(1,9*(Population_Num/10));
        newFitness = zeros(1,9*(Population_Num/10));   
        for Anti = 1:9*(Population_Num/10)    
            new_anti(Anti) = {reshape(bf(:,Anti),M,N)};
            newFitness(Anti) = Distance(cover, new_anti{Anti}); 
            bMSLL(Anti) = newFitness(Anti);
        end
        %%%%%%%%%%%%%%免疫种群与新生种群合并%%%%%%%%%%%%%%%%%%%
        f1 = [af,bf];
        MSLL1 = [aMSLL,bMSLL];
        [SortFitness,Index] = sort(MSLL1);
        Sortf = f1(:,Index);
        fprintf('\n迭代次数：%d 次\n',gen);
    %     gen = gen + 1;
        bestSortf = reshape(Sortf(:,1),M,N);
        trace(gen) = Distance(cover,bestSortf);
        gen = gen + 1;
        toc;
    end
    stego = reshape(stego,[M,N]);
    %% 输出优化结果
    Bestf = Sortf(:,1);                 %最优变量
    IA_stego = reshape(Bestf,M,N);
    trace(end);                       %最优值
%     figure,plot(trace)
%     xlabel('迭代次数')
%     ylabel('目标函数值')
%     title('亲和度进化曲线')
end

