function [SAD A_scale indexx index] =rd_compare(A_est,libfocus,banduse)
libfocus;
load library
%focus=1:498;
% libfocus=focus;
sample=size(library,2);
% banduse=[3:103,114:147,168:220];

N=size(A_est,2);
library = library(banduse,:);
M = length(banduse);


%library_n=library./(ones(M,1)*mean(library));

SAD=[];
index=[];
A_scale=[]; ind_prev=[];
for i=1:N
    if ~isempty(libfocus)
        lib_focus=libfocus;
        % %         lib_focus
        lib_ind=library_index(1,:);
        %         focus=find(lib_focus(1)==lib_ind):find(lib_focus(end)==lib_ind)
        for qw=1:length(lib_focus),
            %     lib_focus(qw)
            focus(qw)=find( lib_ind==lib_focus(qw));
        end
    end
    scale=pinv([A_est(:,i) ones(M,1)])*library(:,focus);
    test=A_est(:,i)*scale(1,:)+ones(M,1)*scale(2,:); % why? to fix scaling ambiguity???
    c_test=test-ones(M,1)*mean(test); % mean removed
    %c_test=test;
    c_library=library(:,focus)-ones(M,1)*mean(library(:,focus));  % mean removed library
    %c_library=library(:,focus);
    theta=acos(sum(c_test.*c_library)./(sqrt(sum(c_test.^2)).*sqrt(sum(c_library.^2))))*180/pi;
    %sorted = sort(theta);
    [val ind]=min(theta);
    while sum(ismember(ind,ind_prev))==1,
        theta(ind)=10000;
        [val ind]=min(theta);
    end
    SAD=[SAD val];
    index=[index focus(ind)];
    A_scale=[A_scale test(:,ind)]; % appropriately scaled version of the estimated signature
    %index=[index library_index(:,ind)];
    ind_prev= [ind_prev ind];
end
ind_set=library_index(:,index);
indexx=ind_set(1,:);


