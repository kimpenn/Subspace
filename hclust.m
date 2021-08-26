    function [assignments, allTrees] = hclust(differences, allAffine, assignments, combined, tree, allTrees, count, majority)


    if size(allAffine,1) == 2
        return;

    else

        [maximumVal, index] = maxk(differences,2,2);

        temp=maximumVal(:,2);

        [M,I] = max(temp);

        %Index of most similiar clique

        accRatio=0.95;

        %Merge Clique I and Clique index(I,2)
        data = vertcat(combined(assignments{I},:),combined(assignments{index(I,2)},:));

        x_bar=mean(data);


        [coefficients, score, latent] = pca(data);


        coefficients = coefficients';


        accuracy = zeros(size(coefficients,2), 2);

        for j =1:size(coefficients,1)
            accuracy(j,1) = j;
            num = sum(latent(1:j, 1));
            denom = sum(latent);
            accuracy(j, 2) = num/denom;

            if j>1 && accuracy(j,2) > accRatio && accuracy(j-1,2)<accRatio
                estimatedDimension = j;
            elseif j<=1
                estimatedDimension = 2;
            end
        end

        subspace = coefficients(1:estimatedDimension, :);

        affine = vertcat(subspace, x_bar);

        %Orthogonalize vectors by Gram-Schmidt
        affine = orth(affine');

        affine = affine';

        allAffine{I}=affine;

        assignments{I}=vertcat(assignments{I}, assignments{index(I,2)});

        assignments{index(I,2)} = [];

        assignments=assignments(~cellfun('isempty',assignments));

        for i=1:size(allAffine,1)

            if size(allAffine{i},1)>0 && size(affine,1)>0

                A = vertcat(affine(1:size(affine,1)-1,:)', zeros(1, size(affine,1)-1));

                lastcol = vertcat(affine(size(affine,1),:)'/sqrt(1+norm(affine(size(affine,1),:))^2), 1/sqrt(1+norm(affine(size(affine,1),:))^2));

                Ap= horzcat(A, lastcol);

                B = vertcat(allAffine{i}(1:size(allAffine{i},1)-1,:)', zeros(1, size(allAffine{i},1)-1));

                lastcol = vertcat(allAffine{i}(size(allAffine{i},1),:)'/sqrt(1+norm(allAffine{i}(size(allAffine{i},1),:))^2), 1/sqrt(1+norm(allAffine{i}(size(allAffine{i},1),:))^2));


                Bp = horzcat(B, lastcol);

                s = svd(Ap'*Bp);

                %Grassman Distance
                %s2=abs(acos(s)).^2;

                %differences(j,i) = sqrt(sum(s2));


                s2=abs(acos(s));

                %LSA Paper Dissimilarity (exponential of Chordal)
                differences(I,i) = exp(-sum((sin(s2).^2)));

                differences(i,I) = differences(j,i);
            end

        end
        
         allAffine{index(I,2)} = [];

        allAffine=allAffine(~cellfun('isempty',allAffine));

        differences(index(I,2),:) = [];

        differences(:,index(I,2)) = [];

        tree{I} = vertcat(tree{I}, tree{index(I,2)});

        tree{index(I,2)} = [];
        tree=tree(~cellfun('isempty',tree));
        
        count=count+1;
        
        allTrees{count,1} = tree;

        [assignments, allTrees] = hclust(differences, allAffine, assignments, combined, allTrees{count,1}, allTrees, count, majority);
    end


