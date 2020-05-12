% Repeated Random Sampling
% Ex:
% rtt = [0.1,0.2,0.3] (rate of testing to training)
% rep = 100;

function ind = fun_rrsid(data_len,rtt,rep,ver_rate)
rng('shuffle')    ;
s = rng;

for i0 = 1:length(rtt)
    
    for i1=1:rep
        n         = round(data_len*rtt(i0));
        ID        = randperm(data_len);
        ID_test   = ID(1:n);
        ID_train  = ID(n+1:end);
        num_verif = floor(length(ID_train)*ver_rate);
        ID_verif  = ID(n+1:n + num_verif);
        ID_train  = ID(n + num_verif + 1:end);
        
        if data_len>2^64-1
            ind(i0,1).test(:,i1)  = ID_test';
            ind(i0,1).train(:,i1) = ID_train';
            ind(i0,1).verif(:,i1) = ID_verif';
        elseif data_len > 2^16-1 && data_len <= 2^64-1
            ind(i0,1).test(:,i1)  = int64(ID_test');
            ind(i0,1).train(:,i1) = int64(ID_train');
            ind(i0,1).verif(:,i1) = int64(ID_verif');
        elseif data_len > 2^08-1 && data_len <= 2^16-1
            ind(i0,1).test(:,i1)  = int16(ID_test');
            ind(i0,1).train(:,i1) = int16(ID_train');
            ind(i0,1).verif(:,i1) = int16(ID_verif');
        elseif data_len <=2^08-1
            ind(i0,1).test(:,i1)  = int8(ID_test');
            ind(i0,1).train(:,i1) = int8(ID_train');
            ind(i0,1).verif(:,i1) = int8(ID_verif');
        end
            
    end

end

end