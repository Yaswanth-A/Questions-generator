%% QUESTION 1 (LINEAR ALGEBRA)
    
    %Initialization
    clear all
    % Define the 16 statements
    statements = {
        'Each superset of a linearly dependent set is linearly dependent.' 
        'Each subset of a linearly independent set is linearly independent.' 
        'Union of any two linearly dependent sets is linearly dependent.' 
        'Intersection of any two linearly independent sets is linearly independent.' 
        'Each subset of a linearly dependent set is linearly dependent.' 
        'Each superset of a linearly independent set is linearly independent.'
        'Union of any two linearly independent sets is linearly independent.' 
        'Intersection of any two linearly dependent sets is linearly dependent.'
        'If span(A) ∩ span(B) = {0}, then A ∪ B is linearly independent.' 
        'If v1, . . . , vn are linearly independent, then v1, v2 − v1, . . . , vn − v1 are linearly independent.'
        'If v1, . . . , vn span V, then v1, v2 −v1, . . . , vn −v1 span V.'
        'If v1, . . . , vn are linearly independent, then v1, v2 − v1, . . . , vn − v1 are linearly dependent.'
        'span(span(A)) = span(A).'
        'If A ⊆ B, then span(A) ⊆ span(B).'
        'span(A ∩ B) ⊆ span(A) ∩ span(B).'
        'span(A) ∩span(B) ⊆ span(A ∩ B).'
    };

    for i=1:5
    %Title
    fprintf('\nQ1V%d\n ',i);

        % Randomly select 4 statements
        selected_indices = randperm(16, 4);
    
        % Display the selected statements
        disp('Which of the following statements are correct?')
        options = char('a' + (1:numel(selected_indices)) - 1);
        for i = 1:numel(selected_indices)
            disp([char(options(i)) ') ' statements{selected_indices(i)}])
        end
        disp(['e) None of the above'])
    
        % Generate the answer key and explanations for the selected statements
        answer_key = [1 1 1 1 0 0 0 0 1 1 1 0 1 1 1 0]; % Correctness of each statement
        explanations = {
            'Let [p,q] be linearly dependent set so ap+bq = 0 if [p,q,r] is superset we can write ap+bq+(0)r = 0 and atleast one of the a,b is not zero, so the superset is linearly dependent  '
            'Use the notion of linear independence'
            'Use the notion of linear dependence'
            'Use the notion of linear independence'
            'Take any linearly independent set and add the zero vector to it'
            'Take any linearly independent set and add the zero vector to it'
            'consider Set A: {v1, v2} where v1 = [1, 0] and v2 = [0, 1]Set B: {v3, v4} where v3 = [2, 0] and v4 = [0, 2]'
            'consider Set A: {v1, v2} where v1 = [1, 0] and v2 = [0, 1]Set B: {v2, v3} where v2 = [0, 1] and v3 = [1, 1]'
            'Set A: {v1} where v1 = [1, 0] Set B: {v2} where v2 = [0, 1]'
            'Trivial'
            'It is enough to show that the given set of vectors spans the vectors v1, v2, ..., vn.'
            'Trivial'
            'Follows from the fact that L(s) is the subspace'
            'Trivial'
            'Take an element form L(A ∩ B), write it as linear combination of elements of A ∩ B and observe'
            'Take two different non-zero vectors such that linear span is the same space' 

        };
    
        % Create options strings
        options = char('a' + (1:numel(selected_indices)) - 1);
        options = [options 'e'];
    
        
    
        % Check the correctness of the answer
        correct_options = options(answer_key(selected_indices) == 1);
        incorrect_options = options(answer_key(selected_indices) == 0);
    
        % Display the answer
        
        if isempty(correct_options)
            disp('answer is e) None of the above')
        else
            disp(['The correct options are: ' char(correct_options)])
            
        end
    
        % Explanation of the answer
        disp('Explanation:')
        
       
            for i = 1:numel(selected_indices)
                index = selected_indices(i);
                %statement = statements{index};
                %answer = answer_key(index);
                explanation = explanations{index};
                %([char(options(i)) ') ' statement])
                %disp(['   Correct answer: ' num2str(answer)])
                disp([char(options(i)) ') ' explanation])
            end
       
    end




    %% QUESTION 2 (LINEAR ALGEBRA)

    clear all

 for i=1:5
    %Title
    fprintf('\nQ2V%d\n ',i);

    % Define matrices P1 to P6 as a cell array
    k = randi([1,10]);
    P = cell(1, 6);
    P{1} = [k, 0, 0; 0, k, 0; 0, 0, k];
    P{2} = [k, 0, 0; 0, 0, k; 0, k, 0];
    P{3} = [0, k, 0; k, 0, 0; 0, 0, k];
    P{4} = [0, k, 0; 0, 0, k; k, 0, 0];
    P{5} = [0, 0, k; k, 0, 0; 0, k, 0];
    P{6} = [0, 0, k; 0, k, 0; k, 0, 0];
    
    % Define matrix Q
    % Generate a symmetric matrix of order 3 with elements in range [-100, 100]
    A = randi([-100, 100], 3);  % Generate random integers in the range [-100, 100]
    
    % Make the matrix symmetric
    Q = (A + A')/2;
    
    
    % Compute X
    X = zeros(3);
    
    for k = 1:6
        X = X + P{k} * Q * P{k}';
    end
    
    % Calculate the correct answer for the sum of diagonal elements
    tr = trace(Q);
    correctAnswer = k^2 * tr * 6;
    
    % Generate random permutations of array indices
    indices = randperm(4);
    
    % Assign jumbled and unique values to variables
    options = [k^2 * tr *6, k^2 * tr * 3, k^2 * tr * 12, k^2 * tr];
    % Generate random permutation of array indices
    permutedIndices = randperm(numel(options));
    
    % Create a new array with jumbled elements
    jumbledOptions = options(permutedIndices);
    
    % Determine the correct option index
    correctOptionIndex = find(jumbledOptions == correctAnswer);
    
    % Display the matrices and the question
    disp('Let ')
    disp('P1:');
    disp(P{1});
    disp('P2:');
    disp(P{2});
    disp('P3:');
    disp(P{3});
    disp('P4:');
    disp(P{4});
    disp('P5:');
    disp(P{5});
    disp('P6:');
    disp(P{6});
    disp('Q:');
    disp(Q);
    
    disp('and X = 𝚺PₖQ(Pₖ)^T  Where k runs from 1 to 6, ')
    disp('Then what is the sum of diagonal elements of the matrix X?')
    fprintf('(a) %d (b) %d (c) %d (d) %d \n', jumbledOptions(1), jumbledOptions(2), jumbledOptions(3), jumbledOptions(4));
    
    
    
    % Display the correct option
    disp('Correct Option:')
    fprintf('Option (%c)\n', char(correctOptionIndex + 96));
    
    
 end

 % Explanation of the answer
    disp('Explanation:')
    disp('We need use the fact that Trace(AB) = Trace(BA)');
    disp('Trace(𝚺PₖQ(Pₖ)^T ) = Trace(𝚺Q(Pₖ)^TPₖ and We can see that (Pₖ)^TPₖ = (k^2)I');
    disp('Where I is identity matrix');
    disp('So Trace(𝚺PₖQ(Pₖ)^T ) = (k^2)𝚺Trace(Q)')
    disp('= 6*(k^2)*Trace(Q)');


 %% QUESTION 3 (OPTIMIZATION)

clear all
for i = 1:5
    %Title
    fprintf('\nQ3V%d\n ',i);

    syms x;
    syms y;

    a = randi([2, 5]);
    b = randi([1, 6]);
    c = randi([3, 7]);
    d = randi([2, 8]);
    e = randi([1, 7]);
    f = randi([3, 9]);

    %Create a function in x and y of form: ax^3 + bx^2 + cy^2 + dx + ey + f
    % f(x, y) = (ax-b)^2 + cx^3 + (dy+e)^2 - f
    coefficient_x3 = c;
    coefficient_x2 = a*a; 
    coefficient_x = (-2)*a*b;
    coefficient_y2 = d*d;
    coefficient_y = 2*d*e;
    constant_term = b*b + e*e - f;

    %computing gradient of f
    % Let grad(f) = [A ;B]
    grad_f = [2*a*(a*x-b) + 3*c*x*x; 2*d*(d*y + e)];

    %Finding stationary points, for which grad(f) = 0
    %The points are (x1, y0) and (x2, y0) 
    
    x1 = (-a*a + sqrt(a^4 + 6*a*b*c)) / (3*c);
    x2 = (-a*a - sqrt(a^4 + 6*a*b*c)) / (3*c);
    y0 = -(e/d);

    %Computing Hessian, H(x, y)
    H = @(x) [2*a^2 + 6*c*x 0 ; 0 2*d^2];
    H_disp = [2*a^2 + 6*c*x 0 ; 0 2*d^2];
    
    % Define the 16 statements
    statements_opt1 = {
        '(x1,y0) is a point of minima'
        '(x2,y0) is a saddle point'
        '(x1,y0) is a point of maxima'
        '(x1,y0) is a saddle point'
        '(x2,y0) is a point of minima'
        '(x2,y0) is a point of maxima'
    };
    % Randomly select 4 statements
    selected_indices_opt1 = randperm(6, 4);
    options_opt1 = char('a' + (1:numel(selected_indices_opt1)) - 1);
    answer_key_opt1 = [1 1 0 0 0 0]; % Correctness of each statement
    %Question
    disp('Consider a 2-variable function f(x, y) defined as follows.');disp(' ');
    fprintf('f(x,y) = %dx^3 + %dx^2 + %dy^2 - %dx + %dy + %d\n', coefficient_x3, ...
        coefficient_x2, coefficient_y2, -coefficient_x, coefficient_y, constant_term);
    disp('Which of the following options are correct for the given function f(x, y)?')
    fprintf('given points (x1, y0) = (%f , %f ) and (x2,y0) = (%f , %f) \n', x1, y0, x2, y0)

    %Generate options
    for i = 1:numel(selected_indices_opt1)
        disp([char(options_opt1(i)) ') ' statements_opt1{selected_indices_opt1(i)}])
    end
    disp(['e) None of the above'])

    %correct options
    % (x1, y0) is a point of minima
    % (x2, y0) is a saddle point
    % Check the correctness of the answer
        correct_options_opt1 = options_opt1(answer_key_opt1(selected_indices_opt1) == 1);
       
    
        % Display the answer
        
        if isempty(correct_options_opt1)
            disp('answer is e) None of the above')
        else
            disp(['The correct options are: ' char(correct_options_opt1)])
            
        end

    %Explanation of the answer
    disp(' '); disp('Explanation:');
    disp('First we find the stationary points for f using ∇f = 0');
    disp('The gradient of f is');
    disp(grad_f);
    fprintf('The roots of %dx^2 + %dx - %d = 0 are %f, %f\n', 3*c, 2*a*a, 2*a*b, x1, x2);
    fprintf('The roots of %dy + %d = 0 is %d\n', 2*d*d, 2*d*e, y0);
    fprintf('Therefore stationary points are (%f, %f) and (%f, %f)\n', x1,y0,x2,y0);
    disp('The hessian of f is ');
    disp(H_disp);disp('');
    fprintf('The hessian of f at (%f, %f) is\n', x1, y0);
    disp(H(x1));disp('');
    fprintf('The eigen values are: %d, %d. So f is minimum at (%f, %f)\n', 2*a^2 + 6*c*x1,2*d^2,x1,y0);
    fprintf('\nThe hessian of f at (%f, %f) is\n', x2, y0);
    disp(H(x2));disp('');
    fprintf('The eigen values of f are: %d, %d. So, (%f, %f) is a saddle point\n',2*a^2 + 6*c*x2,2*d^2,x2,y0);
  
end


%% QUESTION 4 (OPTIMIZATION)

clear all

for i=1:5
    %Title
    fprintf('\nQ4V%d\n', i);

    % Create the function
    f = @(x, y) a*x*x + b*y*y;

    %Set inital guess
    x_int = 1;
    y_int = 1;

    %Generating random values for question; more than 500 variants
    %such that answer lies in [1, 4]
    a = randi([1,5]);
    b = a;

    if (a==1)
        alpha = 0.41 + 0.01*randi([-7, 8]); % Learning rate 
        tolerance = 0.075+0.001*randi([-5,5]); % tolerance on ∇f 
    end

    if (a==2)
        alpha = 0.21+0.01*randi([-3, 4]); % Learning rate 
        tolerance = 0.075+0.001*randi([-5,5]); % tolerance on ∇f 
    end

    if (a==3)
        alpha = 0.14+0.01*randi([-2, 3]); % Learning rate 
        tolerance = 0.055+0.001*randi([-5,5]); % tolerance on ∇f 
    end

    if (a==4)
        alpha = 0.10+0.01*randi([-1, 2]); % Learning rate 
        tolerance = 0.055+0.001*randi([-5,5]); % tolerance on ∇f 
    end

    if (a==5)
        alpha = 0.09+0.001*randi([-10, 8]); % Learning rate 
        tolerance = 0.025+0.001*randi([-5,5]); % tolerance on ∇f 
    end

    %computing gradient of f
    grad_f_x = @(x) 2*a*x; %x-component of gradient
    grad_f_y = @(y) 2*b*y; %y-component of gradient
    
    gradf_x = grad_f_x(1); %At the start 
    gradf_y = grad_f_y(1);

    count = 0; %No. of iterations 
    %Start performing the iterations

    x_temp = x_int;
    y_temp = y_int;
    while ((gradf_x > tolerance) && (gradf_y > tolerance))

        x = x_temp - alpha*gradf_x;
        y = y_temp - alpha*gradf_y; %compute points
        
        gradf_x = grad_f_x(x);
        gradf_y = grad_f_y(y);
        
        x_temp = x; %Updating the kth point of iteration
        y_temp = y;

        count = count + 1; % incrementing the count
    end
    
    % Define the options
    option1_opt2 = count; %correct answer
    option2_opt2 = count+1;
    option3_opt2 = count+3;
    option4_opt2 = count-1;
    
    options_opt2 = [option1_opt2, option2_opt2, option3_opt2, option4_opt2];
    
    
    % Generate random permutation of array indices 
    % We use this for jumbling the options
    permutedIndices_opt2 = randperm(numel(options_opt2));
            
    % Ensuring the options are jumbled each time
    jumbledOptions_opt2 = options_opt2(permutedIndices_opt2);
            
    % Determine the correct option index
    
    correctOptionIndex_opt2 = find(jumbledOptions_opt2 == count);
    
    %Question
    disp('');disp('Suppose an ant is gracing on a surface and the surface is approximated by the function');
    fprintf('f(x,y) = %dx^2 + %dy^2\n',a,b);
    fprintf('The ant is moving on the surface to reach an optimum point. Initially, starting from the point (%d,%d).\n',x_int,y_int);
    disp('The movement of ant is in accordance with the steepest decent algorithm learned.');
    disp('Find out the number of times the ant changes the place from one point to another.');
    fprintf('The ant changes the place(point) until the tolerance of %.3f is obtained on the steepest gradient taken by it at that point\n', tolerance);
    fprintf('Use the leraning factor as %.2f\n\n', alpha);

    %Options 
    disp('Options:');
    fprintf('(a) %d\n', jumbledOptions_opt2(1));
    fprintf('(b) %d\n', jumbledOptions_opt2(2));
    fprintf('(c) %d\n', jumbledOptions_opt2(3));
    fprintf('(d) %d\n\n', jumbledOptions_opt2(4));


    %Correct answer - count

    % Display the correct option
    fprintf('Correct Option : Option (%c)\n\n', char(correctOptionIndex_opt2 + 96));

    %Explanation
    disp('Explanation:');disp('');
    disp('The steepest decent algorithm states that x_{k+1} = x_{k} - α.∇f|_{xk,yk}');disp('');
    disp('Using that algorithm we get the iterated points as:');disp('');
    

    %Print the nth iterated values 
    gradf_x = grad_f_x(1); %At the start 
    gradf_y = grad_f_y(1);
    temp = 0;
    x0 = 1; 
    y0 = 1;
    while (temp ~= count)

        x = x0 - alpha*gradf_x;
        y = y0 - alpha*gradf_y; %compute points
        
        gradf_x = grad_f_x(x);
        gradf_y = grad_f_y(y);
        
        x0 = x; %Updating the kth point of iteration
        y0 = y;
        temp = temp + 1; % incrementing the count
        fprintf('%dth iteration gives (%f,%f) and ∇f at that point is (%f,%f)\n',temp,x,y,gradf_x,gradf_y);
    end
    fprintf('We see that at this point we reached the tolerance for ∇f i.e, %.3f\n', tolerance);
    fprintf('Hence, number of iterations is %d\n', count);
end


%% QUESTION 5 (STATISTICS)

clearvars;  format compact;

for i = 1:5
    %Title
    fprintf('\nQ5V%d\n ',i);
    %Values are choosen so as to ensure the question is hand-solvable
    % Total number of balls
    n = randi([5,15]);
    %Getting random integer coefficients for the linear inequality
    k = randi([1,n]);%constant term in the inequality
    p = randi([1,3]);%coefficient of x
    q = randi([p,4]);%coefficient of y
    
    % Initialize the counter for satisfying cases
    count = 0;
      
    % Iterate over all possible values of x and y
    for x = 1:n
        for y = 1:n
            % Check if the inequality holds
            if (p*x- q*y + k) > 0
                % Increment the counter
                count = count + 1;
            end
        end
    end
    
    % Calculate the probability
    probability = count / (n * n);
    %rounding off to 3 decimal places
    probability = round(probability, 3);
    
    % Generate three random values in the range [0, 1]
    % We use these as options
    options_prob1 = rand(1, 3);
    
    % Round off the options to three decimal places
    rounded_options = round(options_prob1, 3);
    rounded_options = [rounded_options probability];
    
    % Generate random permutation of array indices 
    % We use this for jumbling the options
    permutedIndices_prob1 = randperm(numel(rounded_options));
        
    % Ensuring the options are jumbled each time
    jumbledOptions_prob1 = rounded_options(permutedIndices_prob1);
        
    
    % Determine the correct option index
    correctAnswer = probability;
    correctOptionIndex_prob1 = find(jumbledOptions_prob1 == correctAnswer);
    
    %The Question statement
    fprintf('Assume there are %d balls in a bag of same size and color,', n);
    fprintf(' numbered 1,2,....,%d, were put into a packet. ', n);
    disp('Now Arjun draws a ball from the packet, noted that it is of number "x" and puts back it.');
    disp('Then Bunny also draws a ball from the packet and noted that it is of number "y".');
    fprintf('Then the probability for the inequality %d*x - %d*y + %d > 0 to hold', p,q,k);
    disp('(Round off answer to 3 decimal places)');
    fprintf('\n');
    
    %Displaying the options
    fprintf('(a) %.3f (b) %.3f (c) %.3f (d) %.3f \n\n\n', jumbledOptions_prob1(1), jumbledOptions_prob1(2), jumbledOptions_prob1(3), jumbledOptions_prob1(4));
    
    % Display the correct option
    fprintf('Correct Option : Option (%c)\n\n', char(correctOptionIndex_prob1 + 96));
    
    
end

% Explanation of the answer
    disp('Explanation:')
    disp('Let us take an example problem and solve...All can be solved in the similar way');
    disp(' Let us take the case of x-2y+10 > 0 and there are 9 balls ');
    disp('The total number of possible events is 9*9 = 81');
    disp('From x-2y+10 > 0 we get 2y < x+10.');
    disp('We find that when y=1,2,3,4,5, x can take any value in 1,2,..,9 to make the ineuality hold. Then we have 9*5 = 45 admissible events.' );
    disp('When y = 6, x can be 3,4,..,9 and there are 7 admissible events.');
    disp('When y = 7, x can be 5,6..,9 and there are 5 admissible events.');
    disp('When y = 8, x can be 7,8,9 and there are 3 admissble events.');
    disp('When y = 9, x can be 9 and there is 1 admissible event.');
    disp('So the required probability is (45+7+5+3+1)/81 = 61/81');
    disp('Any question of similar model can be done in the similar way');

%% QUESTION 6 (STATISTICS)

clearvars;  format compact;
for i = 1:5
    %Title
    fprintf('\nQ6V%d\n ',i);
    % Define the winnings and losses
    m = randi([1, 10]); % Winnings for two tails
    n = randi([1, 10]); % Winnings for two heads
    p = randi([1, 10]); % Loss for one head and one tail
    
    % Generate random values for a and b between 0 and 1
    a = rand();
    b = rand();
    
    % Calculate the expected winnings
    E = (1 - a) * (1 - b) * m + a * b * n + (1 - a) * b * (-p) + a * (1 - b) * (-p);
    
    % Calculate the variance
    V = (1 - a) * (1 - b) * (m - E)^2 + a * b * (n - E)^2 + (1 - a) * b * (-p - E)^2 + a * (1 - b) * (-p - E)^2;
    
    % Calculate the ratio V/E
    ratio = V / E;
    
    % Define the options
    option1 = ratio;
    option2 = ratio*2 + ratio/2;
    option3 = ratio*2- ratio/2;
    option4 = ratio*0.65;
    
    options_prob2 = [option1, option2, option3, option4];
    options_prob2 = round(options_prob2, 4);
    
    % Generate random permutation of array indices 
    % We use this for jumbling the options
    permutedIndices_prob2 = randperm(numel(options_prob2));
            
    % Ensuring the options are jumbled each time
    jumbledOptions_prob2 = options_prob2(permutedIndices_prob2);
            
    % Determine the correct option index
    ratio = round(ratio, 4);
    correctOptionIndex_prob2 = find(jumbledOptions_prob2 == ratio);
    
    
    % Display the question and options
    fprintf('You and a friend play a game where you each toss a coin. The coin is biased. The probability of getting a head on the first coin is %f and on the second coin is %f.\n', a, b);
    fprintf('If the upper faces on the coins are both tails, you win $%d; if the faces are both heads, you win $%d; if one shows a head and the other a tail, you lose $%d.\n', m, n, p);
    disp('Your expected winnings is E and variance is V');
    disp('Then what is the V/E (rounded off to 4 decimal places)');
    disp('Options:');
    fprintf('(a) %.4f\n', jumbledOptions_prob2(1));
    fprintf('(b) %.4f\n', jumbledOptions_prob2(2));
    fprintf('(c) %.4f\n', jumbledOptions_prob2(3));
    fprintf('(d) %.4f\n', jumbledOptions_prob2(4));
    
    
    
    % Display the correct option
    fprintf('Correct Option : Option (%c)\n\n', char(correctOptionIndex_prob2 + 96));

   
end
 
 % Explanation of the answer
    disp('Explanation:')
    disp('If the probabilities of getting a head for one coin are denoted as a and for the other coin as b, we can adjust the calculations for the expected winnings and variance accordingly.The updated expected winnings (E) can be calculated as:E = (probability of two tails * winnings for two tails) + (probability of two heads * winnings for two heads) + (probability of one head and one tail * loss)');
    disp('E = ((1 - a) * (1 - b) * m) + (a * b * n) + ((1 - a) * b * (-p)) + (a * (1 - b) * (-p));')
    disp('V = ((1 - a) * (1 - b) * (m - E)^2) + (a * b * (n - E)^2) + ((1 - a) * b * (-p - E)^2) + (a * (1 - b) * (-p - E)^2);')
    disp('m is winnings for two tails, n is winnings for two heads p is loss for one head and one tail');