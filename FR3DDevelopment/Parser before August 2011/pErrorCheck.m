<div class="moz-text-flowed" style="font-family: -moz-fixed">%pErrorCheck.m %covers some of the errors that a user may make when creating the model for a new motif
%Started by Ali Mokdad

for n=1:n %the previously saved value for n is the last node
    if strcmp(Node(n).type,'Motif')
        m=n;
    end
end

if min(size(Node(m+1).P)) ~=1
    fprintf('%s\n','The node after the motif must have: Node(n).P = ones(1,1)*Score;')
end
if sum(Node(m).LIP)~=1
    fprintf('%s\n','Motif error: Node(n).LIP does not add up to 1')
end
if sum(Node(m).RIP)~=1
    fprintf('%s\n','Motif error: Node(n).RIP does not add up to 1')
end
if length(Node(m).Left(:,1))~=numel(Node(m).LIP)
    fprintf('%s\n','Motif error: Node(n).LIP must have equal elements as possible Left bases')
end
if length(Node(m).Right(:,1))~=numel(Node(m).RIP)
    fprintf('%s\n','Motif error: Node(n).RIP must have equal elements as possible Right bases')
end
if length(Node(m).Score(1,1,:))~=length(Node(m).IBases(:,1))
    fprintf('%s\n','Motif error: Number of elements in Node(n).IBases must be the same as Node(n).Score')
end
</div>