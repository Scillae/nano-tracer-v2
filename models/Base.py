# import numpy as np
class Base:
    def __init__(self, base_id, base_type, prev_id, next_id, position,
                 backbone, normal, velocity, angular_velocity, strand_id=-1):
        """
        Base refer to a nucleotide
        :param base_id: id of this base (in a strand). '-1' if not exist
        :param base_type: type of the nitrogen-containing base part of the nucleotide
        :param prev_id: id of the previous base (in a strand). '-1' if not exist
        :param next_id: id of the next base (in a strand). '-1' if not exist
        :param position: where the base locates at
        :param backbone: backbone-to-base unit vector
        :param normal: normal unit vector
        :param velocity: velocity vector
        :param angular_velocity: angular velocity vector
        :param strand_id: id of the strand this base belongs to
        """
        self.base_id = base_id
        self.base_type = base_type
        self.prev_id = prev_id
        self.next_id = next_id
        self.position = position
        self.backbone = backbone
        self.normal = normal
        self.velocity = velocity
        self.angular_velocity = angular_velocity
        self.strand_id = strand_id

    def set_position(self, position):
        """
        set position of base
        :param position: where the base locates at
        :return:
        """
        assert len(position) == 3
        self.position = position
        return self.position

    @staticmethod
    def parse_string(s):
        pass

    @staticmethod
    def parse_list(params):
        """
        unpack list/tuple to Base instance
        :param params: list/tuple of parameters
        :return: Base instance
        """
        (
            position,
            backbone,
            normal,
            velocity,
            angular_velocity,
            base_id,
            strand_id,
            base_type,
            prev_id,
            next_id,
        ) = params
        return Base(base_id, base_type, prev_id, next_id, position,
                    backbone, normal, velocity, angular_velocity, strand_id)
