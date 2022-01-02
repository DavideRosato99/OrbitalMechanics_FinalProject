function [state] = gaplotsafe(~,state,flag)

fid = fopen('gaprogress.txt','a+');
for i = 1:length(state.Score)
   fprintf(fid,'%d %d %d %5.10f\n',state.Population(i,1),state.Population(i,2),state.Population(i,3),state.Score(i)); 
   
end
fclose(fid);



































