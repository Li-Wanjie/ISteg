function [IA_stego,trace] = IA_Post(cover_image,stego_image)
%IA_POST �����㷨����дͼ��������ߴ������
%   ����Ӧ�õ���д����֮��Ϻõ�ͼ��

    tic;
    cover = double(cover_image);
    stego = double(stego_image);
    [M, N] = size(stego);
    imgRes = stego - cover;                     % ���������֮��ĵĲв�
    stego = reshape(stego, [M*N,1]);
    imgRes = reshape(imgRes,[M*N,1]);
    modifyRange = find(imgRes == 1 | imgRes == -1);      % ȡ����ͼ���޸ĵ������ֵ

    Anti_D = M * N;                             % ���߸���ά��,M��N��ͼ��ߴ��С
    Population_Num = 30;                       % ���߸�����Ŀ(��Ⱥ�п�������)
    Anti_Gengs_Upper = 255;                     % ��������ȡֵ���Ͻ�
    Anti_Gengs_Lower = 0;                       % ��������ȡֵ���½�
    IA_Maximum = 100;                           % ������ߴ���,һ��ȡֵ��Χ100-500
    gen = 1;                                    % ���ߴ���
    Ncl = 40;                                   % ��¡����
    C = 4;                                      %  ��д����Ԫ�Ӽ�4����Ӱ��������Ϣ��ȡ
    randNumtemp = 0;

    %% ��ʼ����Ⱥ
    f = repmat(stego, 1, Population_Num);
    cond_anti = cell(1,Population_Num);
    Fitness = zeros(1,Population_Num);

    %% ���㿹��-��ԭ�׺Ͷ�
    for Anti = 1:Population_Num    
        cond_anti(Anti) = {reshape(f(:,Anti),M,N)};
        Fitness(Anti) = Distance(cover, cond_anti{Anti});    
    end
    % ����-��ԭ�׺ͶȰ���������
    [SortFitness, Index]=sort(Fitness);
    Sortf=f(:, Index);

    %% ����ѭ��
    while gen <= IA_Maximum
        for i = 1:Population_Num/10
            %%%%%ѡ����-��ԭ�׺Ͷ�ǰPopulation_Num/10����������߲���%%%%%%
            antibody = Sortf(:,i);
            clone_antibody = repmat(antibody,1,Ncl);        % �����¡��Ncl=10
            for j = 1:Ncl 
                %%%%%%%%%%%%%%%%%����%%%%%%%%%%%%%%%%%%%
                randNum(gen) = modifyRange(randi(numel(modifyRange),1,1));      % ȡ���������ʾ������쿹���1����Ԫ
                if randNum(gen) ~= randNumtemp  % �ж��Ƿ�����һ���������ص�λ����ͬ
                    alfa = randsrc(1,1,[-1 1]); % ����-1��1�������Ƕ��ĵ�Ԫ���мӼ�
                    if (clone_antibody(randNum(gen),j) + alfa * C) >= Anti_Gengs_Lower & (clone_antibody(randNum(gen),j) + alfa * C) <= Anti_Gengs_Upper
                        clone_antibody(randNum(gen),j) = clone_antibody(randNum(gen),j) + alfa * C;
                    end
                end
                randNumtemp = randNum(gen);
            end
            clone_antibody(:,1) = Sortf(:,i);               %������¡Դ����
            %%%%%%%%%%��¡���ƣ������׺Ͷ���ߵĸ���%%%%%%%%%%
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
        %%%%%%%%%%%%%%%%%%%%%%%��Ⱥˢ��%%%%%%%%%%%%%%%%%%%%%%%
        bf = repmat(stego, 1, 9*(Population_Num/10));   % �¿������
        for w = 1:9*(Population_Num/10) 
            %%%%%%%%%%%%%%%%%�������%%%%%%%%%%%%%%%%%%%
            randNum1 = modifyRange(randperm(numel(modifyRange),gen));      % ȡ���������ʾ������쿹���gen����Ԫ
            for q = 1:gen
                alfa1 = randsrc(1,1,[-1 1]);           % ����-1��1�������Ƕ��ĵ�Ԫ���мӼ�
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
        %%%%%%%%%%%%%%������Ⱥ��������Ⱥ�ϲ�%%%%%%%%%%%%%%%%%%%
        f1 = [af,bf];
        MSLL1 = [aMSLL,bMSLL];
        [SortFitness,Index] = sort(MSLL1);
        Sortf = f1(:,Index);
        fprintf('\n����������%d ��\n',gen);
    %     gen = gen + 1;
        bestSortf = reshape(Sortf(:,1),M,N);
        trace(gen) = Distance(cover,bestSortf);
        gen = gen + 1;
        toc;
    end
    stego = reshape(stego,[M,N]);
    %% ����Ż����
    Bestf = Sortf(:,1);                 %���ű���
    IA_stego = reshape(Bestf,M,N);
    trace(end);                       %����ֵ
%     figure,plot(trace)
%     xlabel('��������')
%     ylabel('Ŀ�꺯��ֵ')
%     title('�׺ͶȽ�������')
end

