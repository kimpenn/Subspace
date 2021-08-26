    % mat dataset = matrix input (row: cell, column: gene)
    % csvfile labels = true labels of samples given by dataset authors
    % int k = number of clusters
    % csvfile edges = strength output of SNN-Cliq from SNN.m (https://github.com/BIOINSu/SNN-Cliq/blob/master/source%20code/SNN.m)
    % csvfile cliques = cluster output of SNN-Cliq from Cliq.py (https://github.com/BIOINSu/SNN-Cliq/blob/master/source%20code/Cliq.py)

    function [groups]=subspaceClustering(dataset, labels, numSubspaces, edgeFile, cliqueFile)

    %===Process data and labels===%
    allCellNames = unique(labels);
    allData = cell(size(allCellNames, 1));

    for i=1:size(allCellNames,1)
        x=find(labels == string(allCellNames{i}));
        allData{i} = dataset(x,:);
    end

    endpoints = zeros(numSubspaces,2);

    for i=1:numSubspaces-1
        endpoints(1,1) = 1;
        endpoints(1,2) = size(allData{1},1);
        endpoints(i+1, 1) = endpoints(i,2) + 1;
        endpoints(i+1, 2) = endpoints(i,2) + size(allData{i+1},1);
    end

    combined = [];
    for i=1:numSubspaces
        combined = vertcat(combined, allData{i});
    end

    %===Process output of SNN-Cliq===%
    edges = importdata(edgeFile);
    cliques=importdata(cliqueFile);
    numCliques = 1000000;
    assignments = cell(numCliques,1);

    for j=1:numCliques
        assignments{j}=zeros(size(cliques));
    end

    for i=1:size(cliques,1)
        if cliques(i,1)>0
            assignments{cliques(i,1)}(i,1)=i;
        end
    end

    for j=1:numCliques
        assignments{j}(assignments{j}==0)=[];
    end

    assignments=assignments(~cellfun('isempty',assignments));
    numCliques = size(assignments,1);
    majority = zeros(size(assignments,1),1);

    for i=1:size(assignments,1)
        tallies=zeros(numSubspaces,1);
        for j=1:size(assignments{i},1)
            for k=1:size(endpoints,1)
                if assignments{i}(j,1)>=endpoints(k,1) && assignments{i}(j,1)<=endpoints(k,2)
                    tallies(k,1) = tallies(k,1)+1;
                end

                if tallies(k,1)/size(assignments{i},1) >0.6
                    majority(i,1) = k;
                end
            end
        end
    end

    %===Compute subspaces for cliques===%
    estimatedDimensions = zeros(numCliques,1);
    estimatedDimensionsSVD = zeros(numCliques,1);
    sampleMeans = cell(numCliques,1);
    allAffine = cell(numCliques,1);
    accRatio=0.95;

    for i=1:numCliques

        data=combined(assignments{i},:);

        x_bar=mean(data);
        sampleMeans{i}=x_bar;

        [r, col] = size(data);

        newdata = zeros(r, col);

        for j=1:r
            newdata(j, :) = data(j, :) - x_bar;
        end

        if size(data,1)>2

            result=svds(data,size(data,1)-1);

            [res_x, idx_of_result] = knee_pt(result,1:size(data,1)-1,1);

            estimatedDimensionsSVD(i,1) = res_x;

            [coefficients, score, latent] = pca(data);


            coefficients = coefficients';

            accuracy = zeros(size(coefficients,2), 2);

            for j =1:size(coefficients,1)
                accuracy(j,1) = j;
                num = sum(latent(1:j, 1));
                denom = sum(latent);
                accuracy(j, 2) = num/denom;

                if j>1 && accuracy(j,2) > accRatio && accuracy(j-1,2)<accRatio
                    estimatedDimensions(i,1) = j;
                elseif j<=1
                    estimatedDimensions(i,1) = 2;
                end
            end

            allRecAcc{i} = accuracy;

            allCoefficients{i} = coefficients;
            allSubspaces{i} = coefficients(1:estimatedDimensions(i,1), :);
            affine = vertcat(allSubspaces{i}, sampleMeans{i});

            %Orthogonalize vectors by Gram-Schmidt
            affine = orth(affine');
            allAffine{i}=affine';
        end
    end

    %===Create the affinity matrices===%
    %Compute smallest principal angles between cliques
    differences = zeros(numCliques, numCliques);
    allAngles = cell(numCliques, numCliques);
    allSmallest = zeros(numCliques, numCliques);
    for j=1:numCliques
        for i=j:numCliques
            if size(allAffine{i},1)>0 && size(allAffine{j},1)>0

                A = vertcat(allAffine{j}(1:size(allAffine{j},1)-1,:)', zeros(1, size(allAffine{j},1)-1));

                lastcol = vertcat(allAffine{j}(size(allAffine{j},1),:)'/sqrt(1+norm(allAffine{j}(size(allAffine{j},1),:))^2), 1/sqrt(1+norm(allAffine{j}(size(allAffine{j},1),:))^2));

                Ap= horzcat(A, lastcol);

                B = vertcat(allAffine{i}(1:size(allAffine{i},1)-1,:)', zeros(1, size(allAffine{i},1)-1));

                lastcol = vertcat(allAffine{i}(size(allAffine{i},1),:)'/sqrt(1+norm(allAffine{i}(size(allAffine{i},1),:))^2), 1/sqrt(1+norm(allAffine{i}(size(allAffine{i},1),:))^2));

                Bp = horzcat(B, lastcol);

                s = svd(Ap'*Bp);

                allAngles{j,i} = s;
                allSmallest(j,i) = radtodeg(abs(acos(s(1,1))));
                allSmallest(i,j) = allSmallest(j,i);


                s2=abs(acos(s));

                %LSA Paper Dissimilarity (exponential of Chordal)
                differences(j,i) = exp(-sum((sin(s2).^2)));


                differences(i,j) = differences(j,i);
            end
        end
    end

    %===Compute nearest neighbors for intersection matrix===%
    distances = zeros(size(combined,1), numCliques);

    for i=1:numCliques
        shiftedPoints = zeros(size(combined));
        projectedPoints = zeros(size(combined));
        closest = zeros(size(combined));
        for j=1:size(combined,1)
            shiftedPoints(j,:) = combined(j,:) - sampleMeans{i};
            for k=1:size(allSubspaces{i},1)
                closest(j,:) = closest(j,:)+(dot(shiftedPoints(j,:),allSubspaces{i}(k,:))/dot(allSubspaces{i}(k,:),allSubspaces{i}(k,:)))*allSubspaces{i}(k,:);
            end
            closest(j,:) = closest(j,:) + sampleMeans{i};

            comp = vertcat(combined(j,:), closest(j,:));
            distances(j,i) = pdist(comp,'cityblock');
        end
    end

    %Cutoffs with cusum
    %https://www.mathworks.com/help/signal/ref/cusum.html
    allCuSum= zeros(numCliques,1);
    cliqueNeighbors = cell(numCliques, 1);
    for i=1:numCliques
        [B,indexes] = sort(distances(:,i));
        ipx = cusum(B(:,1));
        allCuSum(i,1) = ipx;
        %cutoffs(i,1) = res_x;
        [m,indices] = mink(distances(:,i), ipx);
        cliqueNeighbors{i} = indices;
    end


    intersectionMat=zeros(numCliques, numCliques);
    for i=1:numCliques
        for j=1:numCliques
            intersectionMat(i,j) = length(intersect(cliqueNeighbors{i},cliqueNeighbors{j}));

            intersectionMat(j,i) = intersectionMat(i,j);
        end
    end

    %===Normalize and combined affinity matrices===%
    intersectionMatTemp=intersectionMat/max(intersectionMat(:));
    combinedDistMatrix=(differences+2*intersectionMatTemp)/3;

    %===Perform spectral clustering===%
    numSubspaces = size(allData, 1);
    group=SpectralClustering(combinedDistMatrix, numSubspaces);

   
    trueClusters100=cell(numSubspaces,1);
    trueClustersIndex100 = cell(numSubspaces,1);
    for i=1:numSubspaces
        x=find(group==i);
        trueClusters100{i} = majority(x,:);
        trueClustersIndex100{i} = x;
    end

    numClusters = numSubspaces;
    input=combined;

    %===Merge cliques and re-compute subspaces===%
    allData=cell(numClusters,1);
    for i=1:numClusters
        allData{i} = [];
        for j=1:size(trueClustersIndex100{i},1)
            allData{i} = vertcat(allData{i},input(assignments{trueClustersIndex100{i}(j,1)},:));
        end
    end


    accRatio=0.95;
    allRecAcc = cell(numClusters,1);
    allSubspaces=cell(numClusters,1);
    allPoints=cell(numClusters,1);
    allCoefficients=cell(numClusters,1);
    estimatedDimensions = zeros(numClusters,1);
    estimatedDimensionsSVD = zeros(numClusters,1);
    sampleMeans = cell(numClusters,1);
    allAffine = cell(numClusters,1);


    for i=1:numClusters

        data=allData{i};

        x_bar=mean(data);
        sampleMeans{i}=x_bar;

        [r, col] = size(data);

        newdata = zeros(r, col);

        for j=1:r
            newdata(j, :) = data(j, :) - x_bar;
        end

        if size(data,1)>2

            result=svds(data,size(data,1)-1);

            [res_x, idx_of_result] = knee_pt(result,1:size(data,1)-1,1);

            estimatedDimensionsSVD(i,1) = res_x;

            [coefficients, score, latent] = pca(data);


            coefficients = coefficients';

            accuracy = zeros(size(coefficients,2), 2);

            for j =1:size(coefficients,1)
                accuracy(j,1) = j;
                num = sum(latent(1:j, 1));
                denom = sum(latent);
                accuracy(j, 2) = num/denom;

                if j>1 && accuracy(j,2) > accRatio && accuracy(j-1,2)<accRatio
                    estimatedDimensions(i,1) = j;
                elseif j<=1
                    estimatedDimensions(i,1) = 2;
                end
            end

            allRecAcc{i} = accuracy;

            allCoefficients{i} = coefficients;

            allSubspaces{i} = coefficients(1:estimatedDimensions(i,1), :);

            affine = vertcat(allSubspaces{i}, sampleMeans{i});

            %Orthogonalize vectors by Gram-Schmidt
            affine = orth(affine');
            allAffine{i}=affine';
        end
    end

    distances = zeros(size(input,1), numClusters);


    for i=1:numClusters
        shiftedPoints = zeros(size(input));
        projectedPoints = zeros(size(input));
        closest = zeros(size(input));
        for j=1:size(input,1)

            shiftedPoints(j,:) = input(j,:) - sampleMeans{i};
            for k=1:size(allSubspaces{i},1)
                closest(j,:) = closest(j,:)+(dot(shiftedPoints(j,:),allSubspaces{i}(k,:))/dot(allSubspaces{i}(k,:),allSubspaces{i}(k,:)))*allSubspaces{i}(k,:);
            end
            closest(j,:) = closest(j,:) + sampleMeans{i};

            comp = vertcat(input(j,:), closest(j,:));
            distances(j,i) = pdist(comp,'cityblock');
        end
    end


    [M,I] = min(distances');
    groups=I';

    return
