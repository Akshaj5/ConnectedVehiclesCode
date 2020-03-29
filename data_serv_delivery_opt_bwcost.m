function [bw_cost,mem] = data_serv_delivery_opt_bwcost(A0,A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14,A15,A16,A17,A18,A19,A20,A21,A22,A23,A24,A25,A26,A27,A28,A29,A30,A31,A32,A33,A34,A35,A36,A37,A38,A39,A40,A41,A42,A43,A44,A45,A46,A47,A48,N,M,beta,mem_edge,mem_occup,mem_app,x,v2e_comtime,v2e_trvtime,vel_free,bandwidth,density_jam,density,bw_const,l_cov)

nov = N;
noe = M;

cvx_begin
    variable mem(nov,noe)
    
    expression mem_int(noe);
    expression bytes_received_min(nov,noe);
    expression bw_util(noe);
    expression utility(nov);
    expression utility_sum;
    expression bw_cost;
    expression mem_accum_vehicle(nov);
    expression mem_accum_edge(noe);
    
    for j = 1:noe
        for i = 1:nov
            mem_int(j) = mem_int(j) + mem(i,j);
            bytes_received_min(i,j) = x(i,j) * bandwidth(j)/(density_jam(j)*(vel_free(j)/3600)*(1-(density(j)/density_jam(j))));
        end
    end
    
    for j = 1:noe
        bw_util(j) = (((mem_int(j)*(vel_free(j)/3600)*(1-(density(j)/density_jam(j))))/l_cov(j))+bw_const(j))/bandwidth(j);
        bw_cost = bw_cost + beta * (1+bw_util(j))^2;
    end
    
    for i = 1:nov
        for j = 1:noe
            mem_accum_vehicle(i) = mem_accum_vehicle(i) + mem(i,j)*x(i,j);
        end
    end

    for j = 1:noe
        for i = 1:nov
            mem_accum_edge(j) = mem_accum_edge(j) + mem(i,j)*x(i,j);
        end
    end
    
    minimize bw_cost
    subject to
        for i = 1:nov
            mem_accum_vehicle(i) == mem_app(i);
        end
        for i = 1:nov
            for j = 1:noe
                if (x(i,j) == 1)
                    mem(i,j) >= 0;
                end
            end
        end
        for i = 1:nov
            for j = 1:noe
                if (x(i,j) == 1)
                    mem(i,j) <= mem_app(i);
                end
            end
        end
        for i = 1:nov
            for j = 1:noe
                if (x(i,j) == 0)
                    mem(i,j) == 0;
                end
            end
        end
        for j = 1:noe
            mem_accum_edge(j) <= mem_edge(j) - mem_occup(j);
        end
        for i = 1:nov
            for j = 1:noe
                mem(i,j)*v2e_comtime(i,j) <= mem(i,j)*v2e_trvtime(i,j);
            end
        end
        for i = 1:nov
            for j = 1:noe
                bytes_received_min(i,j) >= mem(i,j);
            end
        end
cvx_end