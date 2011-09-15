% Truncate the model

        n = n + 1;
        Node(n).type = 'Hairpin';
        Node(n).MiddleIndex = [a:a];
        Node(n).LeftLetter  = '*';
        Node(n).RightLetter = '';
        Node(n).LeftIndex   = a;
        Node(n).RightIndex  = a-1;

        Node(n).Comment = [ ' // Hairpin node type *' ];

        if Verbose > 0,
          fprintf('%3d Hairpin type *\n', n);
        end

        StarScore = zeros(5,5);
        StarScore(5,5) = 1;

        Node(n).IBases(1,1) = 1;
        Node(n).IBases(1,2) = 1;
        Node(n).SubsProb(:,:,1) = StarScore;
        Node(n).InteractionComment{1} = ' // Hairpin for truncation';

        EndLoop = 1;
