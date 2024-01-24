using LibSerialPort
using Revise
includet("LibPiGPIO.jl")
using .LibPiGPIO
using LibSerialPort
const LP = LibPiGPIO
const LS = LibSerialPort
gpioInitialise()

#1.参数预定义
    init_rot = [576 812] #机械臂初始状态时传感器数值
    init_pul = [1500,1500] #初始高频脉冲宽度
    r_rot = 0.017  #绕线盘绕绳半径/m
    rot_range = 270 #舵机角度范围
    num_servo = 2  #舵机个数
#2.预定义通讯接口
    #serial通讯,用于舵机控制信号传输
    ports = get_port_list()
    serial_handle = LS.open(ports[1], 9600)
     #LS.close(serial_handle)--关闭串口

    #i2c通讯
    i2c_handle = i2cOpen(1,0x14,0)
    #  i2cClose(i2c_handle)

#3.自定义函数
    #---控制单个舵机
        function singal_servo_run(id,tar_angle,time)
            len_data = 8
            sent_data = zeros(UInt8,len_data+2,1)
            sent_data[1:2] = [0x55;0x55]; sent_data[4] = 0x03
            sent_data[3] = len_data #数据长度
            sent_data[5] = 1 #舵机个数
            sent_data[6:7] = [0xff & time; 0xff & (time >> 8)] #时间
            sent_data[8] = id #id
            sent_data[9:10] = [0xff & tar_angle[id]; 0xff & (tar_angle[id] >> 8)]
            write(serial_handle,sent_data)
        end

    #---同时控制多个舵机
        """
        tar_angle为目标舵机脉冲高频宽度，time为舵机转动的时间（即控制舵机转动速度）
        与控制板的通讯协议
        数据包格式
           1    2     3       4     5        6          7          8      9              10   .....
        ||0x55 0x55|数据长度|0x03|舵机个数|time高八位|time低八位| 舵机id|舵机角度高八位|舵机角度低八位|....
        """
        function servo_run(tar_angle,time) #tar_angle::Vector(2) , time--ms
            len_data = num_servo*3 + 5
            sent_data = zeros(UInt8,len_data+2,1)
            sent_data[1:2] = [0x55;0x55]; sent_data[4] = 0x03 #不变
            sent_data[3] = len_data #数据长度
            sent_data[5] = num_servo #舵机个数
            sent_data[6:7] = [0xff & time; 0xff & (time >> 8)] #时间

            for id in 1:num_servo
                sent_data[5+3*id] = id #id
                sent_data[6+3*id:7+3*id] = [0xff & tar_angle[id]; 0xff & (tar_angle[id] >> 8)]
            end
            write(serial_handle,sent_data)
        end
        # servo_run(init_pul,2000) #复位
        # rot_data = zeros(Int16, 1,2)
        # @time i2cReadI2CBlockData(i2c_handle,0x00,rot_data,4)
        # init_rot = rot_data
    #---传感器数值转化为角度
        function sensor_angle(rot_data)
            angle = -(rot_data-init_rot)/1024*360
            return angle
        end
    #----construct_q
        function angle2q(q_angle)
            a = 0.12252
            num_angle = length(q_angle)
            q = zeros(2*num_angle+4,1)
            # q[3:4] = [a[1],0]
            q[3:4] = [a,0]
            for i in 1:num_angle
                q[3+2*i] = q[1+2*i] + a*cos(q_angle[i]) #a[i+1]
                q[4+2*i] = q[2+2*i] + a*sin(q_angle[i])
            end
            q
        end
    #---控制
        """
        del_l为绳长变化量，res_ang为当下舵机脉冲高频宽度
        返回值为新的舵机脉冲高频宽度（即下一步的res_ang）与各关节角度（°）
        """
        function for_lya_control(del_u)
            rot_data = zeros(Int16, 1,length(del_u))
            del_ang = zeros(Int64,num_servo,1)
            for i in 1:length(del_u)
                del_ang[i] = floor(Int,180/pi/r_rot*2000/270*del_u[i])
                if iseven(i)
                    del_ang[i] = -del_ang[i]
                end
            end
            mov_time = 50
            tar_angle = init_pul - del_ang

                if tar_angle[1]>2500
                    tar_angle[1]=2500
                elseif tar_angle[1]<500
                    tar_angle[1]=500
                end
                if tar_angle[2]>2500
                    tar_angle[2]=2500
                elseif tar_angle[2]<500
                    tar_angle[2]=500
                end
            servo_run(tar_angle,mov_time)
            @time LP.time_sleep(mov_time/1000)
            i2cReadI2CBlockData(i2c_handle,0x00,rot_data,2*length(del_u))
            q_angle = sensor_angle(rot_data)
            return q_angle,tar_angle
        end

        
# #4.应用
    # singal_servo_run(1,[500,500],2000)
    #     servo_run([500,500],1000)

        # servo_run(init_pul,2000) #复位
        # rot_data = zeros(Int16, 1,2)
        # @time i2cReadI2CBlockData(i2c_handle,0x00,rot_data,4)
        # rot_data

    #     for i in 1:80
    #         del_l = [0.0002*i,0.0002*i]
    #         tar_angle,rot_angle = for_lya_control(del_l)
    #     end
    
    # starttick = gpioTick()
    # endtick = gpioTick()
    # difftick_mics = convert(Int64, endtick - starttick)
    # difftick_sec = difftick_mics/1000000
    

# angle2q([0,0])