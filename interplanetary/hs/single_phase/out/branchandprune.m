function branchandprune(a,limit)

load 'prune.dat' -ASCII;


counter=1;
for i=1:length(prune)
	if ((prune(i,end-1)==141)&(prune(i,end)==a(1))) 
		prune2(counter,:)=prune(i,:);
		counter=counter+1;
	end
end
N1=counter-1;
for i=1:length(prune)
	if ((prune(i,end-1)==a(1))&(prune(i,end)==a(2))) 
		prune2(counter,:)=prune(i,:);
		counter=counter+1;
	end
end
N2=counter-N1-1;
for i=1:length(prune)
	if ((prune(i,end-1)==a(2))&(prune(i,end)==a(3))) 
		prune2(counter,:)=prune(i,:);
		counter=counter+1;
	end
end
N3=counter-N1-N2-1;
for i=1:length(prune)
	if ((prune(i,end-1)==a(3))&(prune(i,end)==141)) 
		prune2(counter,:)=prune(i,:);
		counter=counter+1;
	end
end
N4=counter-N1-N2-N3-1;
clear prune;
prune=prune2;
clear prune2;

value=0; 
valuemin=0;
disp(' ');
for i=1:N1
        if 1
            disp(sprintf('Remaining: %8.3f %%',(N1-(i-1))/N1*100))
        end
        if (prune(i,1)<5844)
            continue
        end
        if (prune(i,1)>9497)
            continue
        end
        for j=N1+1:N1+N2
            if (prune(i,2)+60 > prune(j,1))
               continue
            end            
            for k=N1+N2+1:N1+N2+N3
                if (prune(j,2)+60 > prune(k,1))
                    continue
                end
                for l=N1+N2+N3+1:N1+N2+N3+N4
                    if (prune(k,2)+60 > prune(l,1))
                        continue
                    end
                    if (prune(l,2) > prune(i,1)+3652.5)
                        continue
                    end
    
                    minim=min([prune(j,1)-prune(i,2),  prune(k,1)-prune(j,2), prune(l,1)-prune(k,2)]);
                    obj=prune(i,3)*prune(j,3)*prune(k,3)*prune(l,3)+minim*0.2/3652.5;
                    if (obj>value) 
                        value=obj;
                        valuemin=minim;
                        if (value>limit)
                            v=[prune(i,1),prune(i,2)-prune(i,1), prune(j,1),prune(j,2)-prune(j,1), prune(k,1),prune(k,2)-prune(k,1),prune(l,1),prune(l,2)-prune(l,1) ];
                            fid= fopen('../../input.dat', 'wt');
                            fprintf(fid, 'Sequence %g %g %g\n', [a(1) a(2) a(3)]);
                            fprintf(fid, '%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n',v);
                            fclose(fid);
                        end
                    end
                end
            end
        end
    end
disp(' ');   
disp(' ');   

if (value<limit)
    disp('Previous solution was better: ');
    type input.dat;
  
else
    disp('Best solution:');
    disp(sprintf('Sequence %g %g %g\n', [a(1) a(2) a(3)]));
    disp(sprintf('%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n',v));
    value
    valuemin
end
